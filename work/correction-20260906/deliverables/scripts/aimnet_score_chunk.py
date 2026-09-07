#!/usr/bin/env python3
"""Smoke-test AIMNet2-2025 interaction-energy scoring on one docked SDF chunk."""

from __future__ import annotations

import argparse
import gzip
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd
import torch
from aimnet_safe_calculator import AIMNet2Calculator
from rdkit import Chem, RDLogger
from scipy.stats import rankdata, spearmanr

RDLogger.DisableLog("rdApp.*")
EV_TO_KCAL_MOL = 23.060547830619


def element_from_atom_name(name: str) -> str:
    letters = "".join(c for c in str(name).upper() if c.isalpha())
    if letters.startswith("CL"):
        return "Cl"
    if letters.startswith("BR"):
        return "Br"
    return letters[:1] or "C"


def load_pocket(path: Path) -> tuple[np.ndarray, np.ndarray]:
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        record = pickle.loads(next(iter(txn.cursor()))[1])
    env.close()
    table = Chem.GetPeriodicTable()
    numbers = np.asarray(
        [table.GetAtomicNumber(element_from_atom_name(x)) for x in record["pocket_atoms"]],
        dtype=np.int64,
    )
    coords = np.asarray(record["pocket_coordinates"], dtype=np.float32)
    return coords, numbers


def load_poses(path: Path) -> list[Chem.Mol]:
    with gzip.open(path, "rb") as handle:
        return [m for m in Chem.ForwardSDMolSupplier(handle, removeHs=False) if m is not None]


def ligand_arrays(mol: Chem.Mol) -> tuple[np.ndarray, np.ndarray, int]:
    conf = mol.GetConformer()
    coords = np.asarray(conf.GetPositions(), dtype=np.float32)
    numbers = np.asarray([a.GetAtomicNum() for a in mol.GetAtoms()], dtype=np.int64)
    charge = int(Chem.GetFormalCharge(mol))
    return coords, numbers, charge


def pack(systems: list[tuple[np.ndarray, np.ndarray, int]]) -> dict[str, np.ndarray]:
    coords, numbers, mol_idx, charges = [], [], [], []
    for i, (xyz, nums, charge) in enumerate(systems):
        coords.append(xyz)
        numbers.append(nums)
        mol_idx.append(np.full(len(nums), i, dtype=np.int64))
        charges.append(charge)
    return {
        "coord": np.concatenate(coords),
        "numbers": np.concatenate(numbers),
        "mol_idx": np.concatenate(mol_idx),
        "charge": np.asarray(charges, dtype=np.float32),
    }


def energy(calc: AIMNet2Calculator, systems: list[tuple[np.ndarray, np.ndarray, int]]) -> np.ndarray:
    with torch.inference_mode():
        result = calc(pack(systems))
    return result["energy"].detach().cpu().numpy().astype(np.float64).reshape(-1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poses", required=True)
    parser.add_argument("--pocket", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()

    pose_path, output_path = Path(args.poses), Path(args.output)
    poses = load_poses(pose_path)
    if not poses:
        raise RuntimeError(f"no readable poses in {pose_path}")
    pocket_coord, pocket_numbers = load_pocket(Path(args.pocket))
    calc = AIMNet2Calculator("aimnet2-2025", device=args.device, compile_model=False)
    pocket_e = float(energy(calc, [(pocket_coord, pocket_numbers, 0)])[0])

    rows: list[dict] = []
    for start in range(0, len(poses), args.batch_size):
        batch = poses[start : start + args.batch_size]
        ligands = [ligand_arrays(m) for m in batch]
        complexes = [
            (
                np.concatenate([pocket_coord, xyz]),
                np.concatenate([pocket_numbers, nums]),
                charge,
            )
            for xyz, nums, charge in ligands
        ]
        ligand_e = energy(calc, ligands)
        complex_e = energy(calc, complexes)
        interaction = (complex_e - pocket_e - ligand_e) * EV_TO_KCAL_MOL
        for mol, le, ce, ie in zip(batch, ligand_e, complex_e, interaction):
            props = mol.GetPropsAsDict()
            mol_id = int(props.get("mol_id", mol.GetProp("_Name").split(":")[-1]))
            smina_affinity = float(props["minimizedAffinity"])
            rows.append(
                {
                    "mol_id": mol_id,
                    "name": mol.GetProp("_Name"),
                    "formal_charge": int(Chem.GetFormalCharge(mol)),
                    "n_atoms": mol.GetNumAtoms(),
                    "ligand_energy_ev": float(le),
                    "complex_energy_ev": float(ce),
                    "pocket_energy_ev": pocket_e,
                    "aimnet_interaction_kcal_mol": float(ie),
                    "aimnet_score": float(-ie),
                    "smina_affinity_kcal_mol": smina_affinity,
                    "smina_score": -smina_affinity,
                }
            )
        print(f"scored {min(start + args.batch_size, len(poses))}/{len(poses)}", flush=True)

    frame = pd.DataFrame(rows)
    library = pd.read_csv(args.library)[["mol_id", "source_index", "label", "paff", "smiles"]]
    frame = frame.merge(library, on="mol_id", how="left", validate="one_to_one")
    frame["aimnet_rank"] = rankdata(-frame.aimnet_score, method="average")
    frame["smina_rank"] = rankdata(-frame.smina_score, method="average")
    frame.sort_values("aimnet_rank").to_csv(output_path, index=False)

    labels = frame.label.to_numpy(int)
    active, decoy = labels == 1, labels == 0
    def auc(scores: np.ndarray) -> float:
        ranks = rankdata(scores, method="average")
        return float((ranks[active].sum() - active.sum() * (active.sum() + 1) / 2) / (active.sum() * decoy.sum()))

    rank_rho = float(spearmanr(frame.aimnet_score, frame.smina_score).statistic)
    summary = {
        "model": "aimnet2-2025",
        "score_definition": "-(E_complex - E_pocket - E_ligand)",
        "units": "kcal/mol",
        "scope": "technical smoke test on one non-random chunk; not full AIMNet2(Score)",
        "n_poses": int(len(frame)),
        "n_actives": int(active.sum()),
        "n_decoys": int(decoy.sum()),
        "pocket_atoms": int(len(pocket_numbers)),
        "aimnet_auroc_within_chunk": auc(frame.aimnet_score.to_numpy()),
        "smina_auroc_within_chunk": auc(frame.smina_score.to_numpy()),
        "aimnet_vs_smina_spearman": rank_rho,
        "aimnet_interaction_min": float(frame.aimnet_interaction_kcal_mol.min()),
        "aimnet_interaction_median": float(frame.aimnet_interaction_kcal_mol.median()),
        "aimnet_interaction_max": float(frame.aimnet_interaction_kcal_mol.max()),
        "output": str(output_path),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("\nTop 10 by AIMNet2 interaction score:")
    print(frame.nsmallest(10, "aimnet_rank")[["mol_id", "label", "paff", "aimnet_interaction_kcal_mol", "smina_affinity_kcal_mol"]].to_string(index=False))


if __name__ == "__main__":
    main()
