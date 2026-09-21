#!/usr/bin/env python3
"""Evaluate the public-checkpoint three-term score reconstruction.

Score = interaction energy + gas-to-CPCM desolvation + local ligand strain.
The ligand is minimized in a fixed protein pocket with the interaction model.
All bound-state terms use that geometry; isolated CPCM minimization supplies
the local-strain reference. CSV interaction values are retained only for comparison. This is not an official
released AIMNet2(Score) implementation. The process-local reference-neighbor
calculator avoids the numerical failures identified in the installed backend.
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from aimnet_safe_calculator import AIMNet2Calculator
from aimnet_interaction_calculator import InteractionCalculator
from aimnet_casf_interaction import read_pocket
from pocket_geometry import check_contacts
from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski
from scipy.stats import rankdata, spearmanr


RDLogger.DisableLog("rdApp.*")
EV_TO_KCAL_MOL = 23.060547830619


def load_poses(path: Path) -> list[Chem.Mol]:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return [m for m in Chem.ForwardSDMolSupplier(handle, removeHs=False) if m]
    return [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m]


def ligand_tensors(mol: Chem.Mol, device: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    coord = torch.as_tensor(
        np.asarray(mol.GetConformer().GetPositions(), dtype=np.float32),
        device=device,
    ).unsqueeze(0)
    numbers = torch.tensor(
        [[atom.GetAtomicNum() for atom in mol.GetAtoms()]],
        dtype=torch.long,
        device=device,
    )
    charge = torch.tensor([Chem.GetFormalCharge(mol)], dtype=torch.long, device=device)
    return coord, numbers, charge


def model_energy(
    model: torch.nn.Module,
    coord: torch.Tensor,
    numbers: torch.Tensor,
    charge: torch.Tensor,
) -> torch.Tensor:
    return model({"coord": coord, "numbers": numbers, "charge": charge})["energy"].reshape(-1)


def gas_energy(
    calc: AIMNet2Calculator,
    coord: torch.Tensor,
    numbers: torch.Tensor,
    charge: torch.Tensor,
) -> float:
    with torch.no_grad():
        result = calc(
            {
                "coord": coord.squeeze(0),
                "numbers": numbers.squeeze(0),
                "charge": charge.float(),
            }
        )
    return float(result["energy"].reshape(-1)[0].detach().cpu())


def fire_minimize(
    model: torch.nn.Module,
    coord: torch.Tensor,
    numbers: torch.Tensor,
    charge: torch.Tensor,
    *,
    fmax: float,
    max_steps: int,
) -> tuple[float, int, float, bool]:
    """Unconstrained single-molecule FIRE minimization in CPCM."""
    x = coord.detach().clone()
    velocity = torch.zeros_like(x)
    dt = 0.1
    dt_max = 0.1
    alpha = 0.1
    n_positive = 0
    max_force = float("inf")

    for step in range(max_steps + 1):
        x.requires_grad_(True)
        energy = model_energy(model, x, numbers, charge)
        force = -torch.autograd.grad(energy.sum(), x)[0]
        max_force = float(force.norm(dim=-1).max().detach().cpu())
        if max_force < fmax:
            return float(energy[0].detach().cpu()), step, max_force, True
        if step == max_steps:
            return float(energy[0].detach().cpu()), step, max_force, False

        with torch.no_grad():
            power = float((velocity * force).sum().cpu())
            if power > 0.0:
                vnorm = velocity.norm().clamp_min(1e-12)
                fnorm = force.norm().clamp_min(1e-12)
                velocity.mul_(1.0 - alpha).add_(alpha * vnorm * force / fnorm)
                n_positive += 1
                if n_positive > 5:
                    dt = min(dt * 1.2, dt_max)
                    alpha *= 0.99
            else:
                velocity.zero_()
                dt *= 0.8
                alpha = 0.1
                n_positive = 0

            velocity.add_(dt * force)
            displacement = dt * velocity
            norm = displacement.norm().clamp_min(1e-12)
            displacement.mul_(min(1.0, 0.1 / float(norm.cpu())))
            x = (x + displacement).detach()

    raise AssertionError("unreachable")


def minimize_in_pocket(calc, coord, numbers, charge, pocket_coord, pocket_numbers,
                       pocket_charge, *, fmax, max_steps):
    """FIRE on ligand coordinates only; return the bound geometry and diagnostics."""
    x = coord.detach().clone()
    velocity = torch.zeros_like(x)
    dt, alpha, n_positive = 0.1, 0.1, 0
    full_numbers = torch.cat((pocket_numbers, numbers.squeeze(0)))
    total_charge = charge.float() + pocket_charge
    initial_energy = None
    for step in range(max_steps + 1):
        result = calc({
            "coord": torch.cat((pocket_coord, x.squeeze(0))).detach().clone(),
            "numbers": full_numbers.clone(), "charge": total_charge.clone(),
        }, forces=True)
        energy = float(result["energy"].reshape(-1)[0].detach().cpu())
        force = result["forces"].reshape(-1, 3)[len(pocket_numbers):].detach().unsqueeze(0)
        if not np.isfinite(energy) or not torch.isfinite(force).all():
            raise ValueError("Nonfinite complex energy or ligand forces")
        if initial_energy is None:
            initial_energy = energy
        max_force = float(force.norm(dim=-1).max().cpu())
        if max_force < fmax or step == max_steps:
            return x, energy, step, max_force, max_force < fmax, initial_energy
        with torch.no_grad():
            power = float((velocity * force).sum().cpu())
            if power > 0:
                velocity.mul_(1-alpha).add_(
                    alpha * velocity.norm().clamp_min(1e-12)
                    * force / force.norm().clamp_min(1e-12))
                n_positive += 1
                if n_positive > 5:
                    dt = min(dt * 1.2, 0.1)
                    alpha *= 0.99
            else:
                velocity.zero_()
                dt *= 0.8
                alpha, n_positive = 0.1, 0
            velocity.add_(dt * force)
            displacement = dt * velocity
            displacement.mul_(min(1.0, 0.1 / float(displacement.norm().clamp_min(1e-12).cpu())))
            x = (x + displacement).detach()
    raise AssertionError("unreachable")


def rotatable_torsions(mol: Chem.Mol) -> list[tuple[int, int, int, int]]:
    """Choose one heavy-atom dihedral for each RDKit rotatable bond."""
    torsions: list[tuple[int, int, int, int]] = []
    for j, k in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
        left = [a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx() != k]
        right = [a.GetIdx() for a in mol.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx() != j]
        if not left or not right:
            continue
        left.sort(key=lambda i: (mol.GetAtomWithIdx(i).GetAtomicNum() == 1, -mol.GetAtomWithIdx(i).GetAtomicNum()))
        right.sort(key=lambda i: (mol.GetAtomWithIdx(i).GetAtomicNum() == 1, -mol.GetAtomWithIdx(i).GetAtomicNum()))
        torsions.append((left[0], j, k, right[0]))
    return torsions


def dihedral_angles(coord: torch.Tensor, torsions: list[tuple[int, int, int, int]]) -> torch.Tensor:
    if not torsions:
        return torch.empty(0, device=coord.device, dtype=coord.dtype)
    idx = torch.tensor(torsions, dtype=torch.long, device=coord.device)
    points = coord[0, idx]
    p0, p1, p2, p3 = points.unbind(dim=1)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / b1.norm(dim=-1, keepdim=True).clamp_min(1e-12)
    v = b0 - (b0 * b1).sum(dim=-1, keepdim=True) * b1
    w = b2 - (b2 * b1).sum(dim=-1, keepdim=True) * b1
    x = (v * w).sum(dim=-1)
    y = (torch.cross(b1, v, dim=-1) * w).sum(dim=-1)
    return torch.atan2(y, x)


def torsion_restrained_refine(
    model: torch.nn.Module,
    mol: Chem.Mol,
    coord: torch.Tensor,
    numbers: torch.Tensor,
    charge: torch.Tensor,
    *,
    force_constant: float,
    max_steps: int,
) -> tuple[torch.Tensor, int, int]:
    """Relax bonds/angles while harmonically preserving rotatable dihedrals."""
    torsions = rotatable_torsions(mol)
    reference = dihedral_angles(coord, torsions).detach()
    x = coord.detach().clone().requires_grad_(True)
    optimizer = torch.optim.LBFGS(
        [x],
        lr=0.5,
        max_iter=max_steps,
        max_eval=max_steps * 2,
        tolerance_grad=2e-3,
        tolerance_change=1e-10,
        history_size=50,
        line_search_fn="strong_wolfe",
    )

    def closure() -> torch.Tensor:
        optimizer.zero_grad(set_to_none=True)
        energy = model_energy(model, x, numbers, charge).sum()
        if torsions:
            delta = dihedral_angles(x, torsions) - reference
            delta = torch.atan2(torch.sin(delta), torch.cos(delta))
            energy = energy + 0.5 * force_constant * delta.square().sum()
        energy.backward()
        return energy

    optimizer.step(closure)
    iterations = int(optimizer.state[x].get("n_iter", 0))
    return x.detach(), iterations, len(torsions)


def molecule_id(mol: Chem.Mol, fallback: int) -> int:
    props = mol.GetPropsAsDict()
    if "mol_id" in props:
        return int(props["mol_id"])
    try:
        return int(mol.GetProp("_Name").split(":")[-1])
    except Exception:
        return fallback


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poses", required=True)
    parser.add_argument("--interaction-csv")
    parser.add_argument("--interaction-kcal", type=float, help="Original interaction value for comparison only")
    parser.add_argument("--pocket", help="Fixed pocket PDB; defaults to pocket.pdb beside poses")
    parser.add_argument("--pocket-charge", type=int, help="Override net pocket charge from PDB")
    parser.add_argument("--complex-max-steps", type=int, default=1000)
    parser.add_argument("--complex-fmax", type=float, default=0.002)
    parser.add_argument("--minimized-poses", help="Output bound poses SDF")
    parser.add_argument("--cpcm-model", required=True)
    parser.add_argument("--gas-model", default="/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt")
    parser.add_argument("--interaction-model", default="aimnet2-2025")
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-poses", type=int, default=250)
    parser.add_argument("--max-steps", type=int, default=1000)
    parser.add_argument("--fmax", type=float, default=0.002)
    parser.add_argument("--refine-steps", type=int, default=300)
    parser.add_argument("--torsion-force-constant", type=float, default=1000.0)
    parser.add_argument("--skip-bound-refinement", action="store_true")
    args = parser.parse_args()
    print("Bound geometry uses fixed-pocket gas FIRE. Legacy CPCM bound-refinement options are ignored.", flush=True)

    poses = load_poses(Path(args.poses))[: args.max_poses]
    ids = [molecule_id(m, i) for i, m in enumerate(poses)]
    if len(ids) != len(set(ids)):
        parser.error("input poses must contain unique mol_id values")
    if args.interaction_csv:
        previous = pd.read_csv(args.interaction_csv)
        if previous["mol_id"].duplicated().any():
            parser.error("interaction CSV must contain one row per mol_id")
        previous_by_id = previous.set_index("mol_id")
        interaction_by_id = previous_by_id["aimnet_interaction_kcal_mol"].to_dict()
    elif args.interaction_kcal is not None:
        previous_by_id = None
        interaction_by_id = None
    else:
        previous_by_id = None
        interaction_by_id = None

    if args.complex_max_steps < 0 or args.complex_fmax <= 0:
        parser.error("complex-max-steps must be nonnegative and complex-fmax positive")
    if args.max_poses <= 0 or not poses:
        parser.error("at least one input pose is required")
    pocket_path = Path(args.pocket) if args.pocket else Path(args.poses).parent / "pocket.pdb"
    pocket_xyz, pocket_z, pdb_charge = read_pocket(pocket_path)
    if len(pocket_z) == 0:
        parser.error("pocket PDB contains no atoms")
    if args.pocket_charge is not None and args.pocket_charge != pdb_charge:
        parser.error("pocket-charge cannot override the validated preparation charge")
    pocket_charge = pdb_charge
    device = args.device
    gas = AIMNet2Calculator(args.gas_model, device=device, compile_model=False)
    interaction = InteractionCalculator(args.interaction_model, device=device, compile_model=False)
    cpcm = torch.jit.load(args.cpcm_model, map_location=device).to(device).eval()
    for parameter in cpcm.parameters():
        parameter.requires_grad_(False)

    output_path = Path(args.output)
    summary_path = Path(args.summary) if args.summary else output_path.with_suffix(".summary.json")
    minimized_path = Path(args.minimized_poses) if args.minimized_poses else output_path.with_suffix(".minimized.sdf")
    for destination in (output_path, summary_path, minimized_path):
        if destination.exists():
            raise FileExistsError(f"Use a new output path; refusing to mix protocols or overwrite: {destination}")
    for destination in (output_path, summary_path, minimized_path):
        destination.parent.mkdir(parents=True, exist_ok=True)
    pocket_coord = torch.as_tensor(pocket_xyz, device=device)
    pocket_numbers = torch.as_tensor(pocket_z, device=device)
    pocket_ev = gas_energy(interaction, pocket_coord.unsqueeze(0), pocket_numbers.unsqueeze(0),
                           torch.tensor([pocket_charge], device=device))
    minimized_writer = Chem.SDWriter(str(minimized_path))
    rows: list[dict] = []
    started = time.perf_counter()
    for index, mol in enumerate(poses):
        mol_id = molecule_id(mol, index)
        row: dict = {
            "mol_id": mol_id,
            "name": mol.GetProp("_Name"),
            "formal_charge": int(Chem.GetFormalCharge(mol)),
            "n_atoms": mol.GetNumAtoms(),
        }
        if previous_by_id is not None and mol_id in previous_by_id.index:
            prior = previous_by_id.loc[mol_id]
            for key in (
                "source_index",
                "label",
                "paff",
                "smiles",
                "ligand_pdb",
                "pose_index",
                "ligunity_rank",
                "ligunity_score",
                "smina_affinity_kcal_mol",
                "smina_score",
            ):
                if key in prior.index:
                    row[key] = prior[key]
        try:
            row["initial_min_heavy_contact_A"] = check_contacts(
                mol.GetConformer().GetPositions(), [a.GetAtomicNum() for a in mol.GetAtoms()],
                pocket_xyz, pocket_z)
            coord, numbers, charge = ligand_tensors(mol, device)
            refined, complex_ev, bound_steps, bound_force, bound_converged, initial_complex_ev = minimize_in_pocket(
                interaction, coord, numbers, charge, pocket_coord, pocket_numbers, pocket_charge,
                fmax=args.complex_fmax, max_steps=args.complex_max_steps)
            refine_steps, n_torsions = 0, len(rotatable_torsions(mol))
            row["final_min_heavy_contact_A"] = check_contacts(
                refined.detach().cpu().numpy(), numbers.detach().cpu().numpy().reshape(-1),
                pocket_xyz, pocket_z)
            with torch.no_grad():
                ecpcm_bound = float(model_energy(cpcm, refined.clone(), numbers, charge)[0].cpu())
            egas = gas_energy(gas, refined.clone(), numbers, charge)
            ecpcm_min, steps, max_force, converged = fire_minimize(
                cpcm,
                refined,
                numbers,
                charge,
                fmax=args.fmax,
                max_steps=args.max_steps,
            )
            initial_interaction_ligand_ev = gas_energy(interaction, coord.clone(), numbers, charge)
            interaction_ligand_ev = gas_energy(interaction, refined.clone(), numbers, charge)
            eint = (complex_ev - pocket_ev - interaction_ligand_ev) * EV_TO_KCAL_MOL
            input_eint = args.interaction_kcal if interaction_by_id is None else interaction_by_id.get(mol_id)
            if not np.isfinite([eint, egas, ecpcm_bound, ecpcm_min]).all():
                raise ValueError("Nonfinite score component")
            desolv = (egas - ecpcm_bound) * EV_TO_KCAL_MOL
            lcse = (ecpcm_bound - ecpcm_min) * EV_TO_KCAL_MOL
            row.update({
                "aimnet_interaction_kcal_mol": eint,
                "input_interaction_kcal_mol": input_eint,
                "pre_minimization_interaction_kcal_mol": (initial_complex_ev - pocket_ev - initial_interaction_ligand_ev) * EV_TO_KCAL_MOL,
                "initial_interaction_ligand_ev": initial_interaction_ligand_ev,
                "interaction_geometry": "fixed_prepared_pocket_aimnet2025_minimized_v3",
                "complex_bound_ev": complex_ev,
                "pocket_ev": pocket_ev,
                "complex_initial_ev": initial_complex_ev,
                "complex_minimization_steps": bound_steps,
                "complex_final_max_force_ev_a": bound_force,
                "complex_optimization_converged": bound_converged,
                "gas_bound_ev": egas,
                "interaction_ligand_ev": interaction_ligand_ev,
                "interaction_model": args.interaction_model,
                "desolvation_gas_model": args.gas_model,
                "cpcm_bound_ev": ecpcm_bound,
                "cpcm_local_min_ev": ecpcm_min,
                "desolvation_kcal_mol": desolv,
                "local_strain_kcal_mol": lcse,
                "aimnet2_composite_kcal_mol": eint + desolv + lcse,
                "optimization_steps": steps,
                "bound_refinement_steps": refine_steps,
                "restrained_rotatable_torsions": n_torsions,
                "final_max_force_ev_a": max_force,
                "optimization_converged": converged,
                "error": "",
            })
            minimized_mol = Chem.Mol(mol)
            for atom_index, position in enumerate(refined.squeeze(0).cpu().numpy()):
                minimized_mol.GetConformer().SetAtomPosition(atom_index, tuple(map(float, position)))
            minimized_mol.SetIntProp("mol_id", mol_id)
            minimized_mol.SetDoubleProp("aimnet_interaction_kcal_mol", eint)
            minimized_mol.SetProp("complex_optimization_converged", str(bound_converged))
            minimized_writer.write(minimized_mol)
            minimized_writer.flush()
            print(
                f"{index + 1}/{len(poses)} mol_id={mol_id} "
                f"Eint={eint:.3f} desolv={desolv:.3f} strain={lcse:.3f} "
                f"score={eint + desolv + lcse:.3f} steps={steps} converged={converged}",
                flush=True,
            )
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
            print(f"{index + 1}/{len(poses)} mol_id={mol_id} ERROR {row['error']}", flush=True)
        rows.append(row)
        pd.DataFrame(rows).to_csv(output_path, index=False)

    minimized_writer.close()
    elapsed = time.perf_counter() - started
    frame = pd.DataFrame(rows)
    for column in ("aimnet2_composite_kcal_mol", "aimnet_interaction_kcal_mol",
                   "desolvation_kcal_mol", "local_strain_kcal_mol",
                   "optimization_converged", "complex_optimization_converged"):
        if column not in frame:
            frame[column] = np.nan
    valid = frame["aimnet2_composite_kcal_mol"].notna()
    frame.loc[valid, "aimnet_interaction_rank"] = rankdata(
        frame.loc[valid, "aimnet_interaction_kcal_mol"], method="average")
    frame.loc[valid, "aimnet2_score_rank"] = rankdata(
        frame.loc[valid, "aimnet2_composite_kcal_mol"], method="average"
    )
    frame.to_csv(output_path, index=False)
    summary = {
        "scope": "reference reconstruction; not an exact reproduction of the unreleased AIMNet2(Score) workflow",
        "score_definition": "Eint + Edesolv + ELCSE; lower is better",
        "interaction_geometry": "fixed_prepared_pocket_aimnet2025_minimized_v3",
        "pocket": str(pocket_path),
        "pocket_charge": pocket_charge,
        "pocket_preparation_record": str(pocket_path.with_suffix(".prep.json")),
        "pocket_preparation": json.loads(pocket_path.with_suffix(".prep.json").read_text()),
        "gas_model": args.gas_model,
        "interaction_model": args.interaction_model,
        "interaction_members": interaction.model_names,
        "cpcm_model": args.cpcm_model,
        "complex_fmax_ev_a": args.complex_fmax,
        "complex_max_steps": args.complex_max_steps,
        "n_complex_optimization_converged": int(frame.loc[valid, "complex_optimization_converged"].fillna(False).sum()),
        "minimized_poses": str(minimized_path),
        "bound_refinement": "AIMNet2(2025) ligand-in-fixed-pocket FIRE; no isolated bound refinement",
        "n_requested": len(poses),
        "n_scored": int(valid.sum()),
        "n_errors": int((~valid).sum()),
        "n_optimization_converged": int(frame.loc[valid, "optimization_converged"].fillna(False).sum()),
        "elapsed_seconds_this_run": elapsed,
        "output": str(output_path),
    }
    if "label" in frame and valid.any():
        evaluation = frame.loc[valid & frame["label"].notna()].copy()
        labels = evaluation["label"].astype(int).to_numpy()
        scores = -evaluation["aimnet2_composite_kcal_mol"].to_numpy(float)
        n_active = int(labels.sum())
        n_total = len(labels)
        n_top = max(1, int(np.ceil(0.01 * n_total)))
        order = np.argsort(-scores, kind="stable")
        active_top = int(labels[order[:n_top]].sum())
        summary.update({
            "n_actives": n_active,
            "n_decoys": n_total - n_active,
            "ef1_percent": float((active_top / n_top) / (n_active / n_total)) if n_active and n_total else None,
            "ef1_top_n": n_top,
            "ef1_actives_in_top_n": active_top,
        })
        if 0 < n_active < n_total:
            ranks = rankdata(scores, method="average")
            n_decoy = n_total - n_active
            summary["auroc"] = float(
                (ranks[labels == 1].sum() - n_active * (n_active + 1) / 2) / (n_active * n_decoy)
            )
        if "smina_score" in evaluation:
            summary["composite_vs_smina_spearman"] = float(
                spearmanr(scores, evaluation["smina_score"].to_numpy(float)).statistic
            )
        summary["composite_vs_interaction_spearman"] = float(
            spearmanr(scores, -evaluation["aimnet_interaction_kcal_mol"].to_numpy(float)).statistic
        )
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {len(frame)} rows to {output_path} in {elapsed:.2f} s")
    print(json.dumps(summary, indent=2))
    print(frame[["mol_id", "aimnet_interaction_kcal_mol", "desolvation_kcal_mol", "local_strain_kcal_mol", "aimnet2_composite_kcal_mol", "optimization_converged"]].to_string(index=False))


if __name__ == "__main__":
    main()
