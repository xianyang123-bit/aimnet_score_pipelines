#!/usr/bin/env python3
"""Calculate AIMNet2 interaction energies for one T3 active-only pose set."""

from __future__ import annotations

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
    xyz, z = [], []
    charge = 0
    table = Chem.GetPeriodicTable()
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        symbol = line[76:78].strip() or "C"
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        z.append(table.GetAtomicNumber(symbol.title()))
        formal_charge = line[78:80].strip()
        if formal_charge:
            charge += int(formal_charge[0]) * (1 if formal_charge[-1] == "+" else -1)
    return np.asarray(xyz, np.float32), np.asarray(z, np.int64), charge


def energies(calc, systems):
    coord = np.concatenate([x for x, _, _ in systems]).astype(np.float32, copy=False)
    numbers = np.concatenate([z for _, z, _ in systems]).astype(np.int64, copy=False)
    mol_idx = np.concatenate([np.full(len(z), i, np.int64) for i, (_, z, _) in enumerate(systems)])
    charge = np.asarray([q for _, _, q in systems], np.float32)
    with torch.inference_mode():
        result = calc({"coord": coord, "numbers": numbers, "charge": charge, "mol_idx": mol_idx})
    return result["energy"].detach().cpu().numpy().astype(np.float64).reshape(-1)


def affinity(mol: Chem.Mol):
    for key in ("minimizedAffinity", "affinity", "vina_affinity"):
        if mol.HasProp(key):
            return float(mol.GetProp(key))
    for key in mol.GetPropNames():
        if "affinity" in key.lower():
            return float(mol.GetProp(key))
    return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pocket", required=True)
    ap.add_argument("--poses", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--batch-size", type=int, default=4)
    args = ap.parse_args()
    poses = [m for m in Chem.SDMolSupplier(args.poses, removeHs=False) if m]
    pocket_xyz, pocket_z, pocket_charge = read_pocket(Path(args.pocket))
    calc = AIMNet2Calculator(args.model, device="cuda", compile_model=False)
    pocket_e = float(energies(calc, [(pocket_xyz, pocket_z, pocket_charge)])[0])
    rows, started = [], time.time()
    for begin in range(0, len(poses), args.batch_size):
        group = poses[begin:begin + args.batch_size]
        systems = []
        for mol in group:
            xyz = np.asarray(mol.GetConformer().GetPositions(), np.float32)
            z = np.asarray([a.GetAtomicNum() for a in mol.GetAtoms()], np.int64)
            ligand_charge = Chem.GetFormalCharge(mol)
            complex_charge = pocket_charge + ligand_charge
            systems.extend([(np.concatenate([pocket_xyz, xyz]), np.concatenate([pocket_z, z]), complex_charge),
                            (xyz, z, ligand_charge)])
        e = energies(calc, systems)
        for i, mol in enumerate(group):
            interaction = (float(e[2*i]) - pocket_e - float(e[2*i+1])) * EV_TO_KCAL_MOL
            rows.append({
                "mol_id": int(mol.GetProp("mol_id")),
                "source_index": int(mol.GetProp("source_index")),
                "layer": mol.GetProp("layer"), "uniprot": mol.GetProp("uniprot"),
                "paff": float(mol.GetProp("paff")), "smiles": mol.GetProp("smiles"),
                "formal_charge": Chem.GetFormalCharge(mol),
                "smina_affinity_kcal_mol": affinity(mol),
                "smina_score": -affinity(mol),
                "aimnet_complex_eV": float(e[2*i]), "aimnet_pocket_eV": pocket_e,
                "aimnet_ligand_eV": float(e[2*i+1]),
                "aimnet_interaction_kcal_mol": interaction,
            })
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(json.dumps({"n_poses": len(poses), "n_scored": len(rows),
                      "pocket_atoms": len(pocket_z), "pocket_formal_charge": pocket_charge,
                      "elapsed_seconds": time.time() - started}, indent=2))


if __name__ == "__main__":
    main()

