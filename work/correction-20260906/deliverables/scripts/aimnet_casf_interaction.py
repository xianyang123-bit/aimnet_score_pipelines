import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from aimnet_safe_calculator import AIMNet2Calculator
from rdkit import Chem, RDLogger

EV_TO_KCAL_MOL = 23.060547830619
RDLogger.DisableLog("rdApp.*")


def read_pocket(path: Path):
    coords, numbers = [], []
    charge = 0
    periodic = Chem.GetPeriodicTable()
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        symbol = line[76:78].strip()
        if not symbol:
            atom_name = "".join(c for c in line[12:16] if c.isalpha()).upper()
            symbol = atom_name[0] if atom_name else "C"
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        numbers.append(periodic.GetAtomicNumber(symbol.title()))
        formal_charge = line[78:80].strip()
        if formal_charge:
            charge += int(formal_charge[0]) * (1 if formal_charge[-1] == "+" else -1)
    return np.asarray(coords, np.float32), np.asarray(numbers, np.int64), charge


def eval_systems(calc, systems):
    coord = np.concatenate([x for x, _, _ in systems]).astype(np.float32, copy=False)
    numbers = np.concatenate([z for _, z, _ in systems]).astype(np.int64, copy=False)
    mol_idx = np.concatenate([np.full(len(z), i, np.int64) for i, (_, z, _) in enumerate(systems)])
    charge = np.asarray([q for _, _, q in systems], np.float32)
    with torch.inference_mode():
        result = calc({"coord": coord, "numbers": numbers, "charge": charge, "mol_idx": mol_idx})
    return result["energy"].detach().cpu().numpy().astype(np.float64).reshape(-1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pocket", required=True)
    parser.add_argument("--poses", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output-all", required=True)
    parser.add_argument("--output-best", required=True)
    parser.add_argument("--output-best-sdf", required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    args = parser.parse_args()

    poses = [m for m in Chem.SDMolSupplier(args.poses, removeHs=False) if m]
    pocket_xyz, pocket_z, pocket_charge = read_pocket(Path(args.pocket))
    print(f"pocket_atoms={len(pocket_z)} poses={len(poses)} cuda={torch.cuda.is_available()}", flush=True)
    calc = AIMNet2Calculator(args.model, device="cuda", compile_model=False)
    pocket_e = float(eval_systems(calc, [(pocket_xyz, pocket_z, pocket_charge)])[0])
    supported = set((calc.metadata or {}).get("implemented_species") or [])
    rows = []
    started = time.time()
    for begin in range(0, len(poses), args.batch_size):
        group = poses[begin:begin + args.batch_size]
        systems, valid = [], []
        for mol in group:
            xyz = np.asarray(mol.GetConformer().GetPositions(), np.float32)
            z = np.asarray([a.GetAtomicNum() for a in mol.GetAtoms()], np.int64)
            unsupported = sorted(set(map(int, z)) - supported) if supported else []
            if unsupported:
                rows.append({"mol_id": mol.GetIntProp("mol_id"), "ligand_pdb": mol.GetProp("ligand_pdb"),
                             "pose_index": mol.GetIntProp("pose_index"), "error": f"unsupported={unsupported}"})
                continue
            ligand_charge = Chem.GetFormalCharge(mol)
            complex_charge = pocket_charge + ligand_charge
            systems.extend([(np.concatenate([pocket_xyz, xyz]), np.concatenate([pocket_z, z]), complex_charge),
                            (xyz, z, ligand_charge)])
            valid.append((mol, ligand_charge))
        if valid:
            energies = eval_systems(calc, systems)
            for i, (mol, q) in enumerate(valid):
                complex_e, ligand_e = map(float, energies[2 * i:2 * i + 2])
                interaction = (complex_e - pocket_e - ligand_e) * EV_TO_KCAL_MOL
                rows.append({
                    "mol_id": mol.GetIntProp("mol_id"), "ligand_pdb": mol.GetProp("ligand_pdb"),
                    "pose_index": mol.GetIntProp("pose_index"), "label": mol.GetIntProp("label"),
                    "ligunity_rank": mol.GetIntProp("ligunity_rank"),
                    "ligunity_score": mol.GetDoubleProp("ligunity_score"),
                    "formal_charge": q, "aimnet_complex_eV": complex_e,
                    "aimnet_pocket_eV": pocket_e, "aimnet_ligand_eV": ligand_e,
                    "aimnet_interaction_kcal_mol": interaction, "error": "",
                })
        if begin == 0 or (begin + len(group)) % 100 == 0:
            print(f"processed={begin + len(group)}/{len(poses)}", flush=True)

    frame = pd.DataFrame(rows)
    frame.to_csv(args.output_all, index=False)
    valid_frame = frame[frame["error"].fillna("") == ""].copy()
    best = valid_frame.loc[valid_frame.groupby("ligand_pdb")["aimnet_interaction_kcal_mol"].idxmin()].copy()
    best = best.sort_values("aimnet_interaction_kcal_mol")
    best["aimnet_interaction_rank"] = np.arange(1, len(best) + 1)
    best.to_csv(args.output_best, index=False)
    selected = set(best["mol_id"].astype(int))
    writer = Chem.SDWriter(args.output_best_sdf)
    for mol in poses:
        if mol.GetIntProp("mol_id") in selected:
            writer.write(mol)
    writer.close()
    summary = {
        "score_definition": "E_complex - E_pocket - E_ligand; lower is better",
        "n_input_poses": len(poses), "n_scored_poses": len(valid_frame),
        "n_ligands": int(frame["ligand_pdb"].nunique()), "n_best_poses": len(best),
        "pocket_atoms": len(pocket_z), "pocket_formal_charge": pocket_charge,
        "elapsed_seconds": time.time() - started,
    }
    Path(args.output_best).with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(best.to_string(index=False))


if __name__ == "__main__":
    main()
