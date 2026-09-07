#!/usr/bin/env python3
"""AIMNet2-2025 interaction-energy smoke test on one completed T3 pose chunk."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from aimnet.calculators import AIMNet2Calculator
from rdkit import Chem, RDLogger
from scipy.stats import rankdata, spearmanr

EV_TO_KCAL_MOL = 23.060547830619
RDLogger.DisableLog("rdApp.*")


def read_pocket(path: Path) -> tuple[np.ndarray, np.ndarray, int]:
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
    return np.asarray(coords, dtype=np.float32), np.asarray(numbers, dtype=np.int64), charge


def read_poses(path: Path, library_path: Path) -> list[dict]:
    library = pd.read_csv(library_path).set_index("mol_id")
    rows = []
    with gzip.open(path, "rb") as handle:
        for mol in Chem.ForwardSDMolSupplier(handle, removeHs=False):
            if mol is None:
                continue
            conf = mol.GetConformer()
            affinity = float(mol.GetProp("minimizedAffinity"))
            name = mol.GetProp("_Name")
            mol_id = int(name.rsplit(":", 1)[-1])
            meta = library.loc[mol_id]
            rows.append({
                "mol_id": mol_id,
                "source_index": int(meta.source_index),
                "label": int(meta.label),
                "paff": float(meta.paff) if not pd.isna(meta.paff) else np.nan,
                "smina_affinity": affinity,
                "charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
                "numbers": np.asarray([atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=np.int64),
                "coords": np.asarray(conf.GetPositions(), dtype=np.float32),
            })
    return rows


def eval_systems(calc: AIMNet2Calculator, systems: list[tuple[np.ndarray, np.ndarray, int]]) -> np.ndarray:
    # Sparse batching avoids zero-padding ghost atoms in the external D3/LR path.
    coord = np.concatenate([xyz for xyz, _, _ in systems], axis=0).astype(np.float32, copy=False)
    numbers = np.concatenate([z for _, z, _ in systems], axis=0).astype(np.int64, copy=False)
    mol_idx = np.concatenate([
        np.full(len(z), i, dtype=np.int64) for i, (_, z, _) in enumerate(systems)
    ])
    charge = np.asarray([q for _, _, q in systems], dtype=np.float32)
    with torch.inference_mode():
        result = calc({"coord": coord, "numbers": numbers, "charge": charge, "mol_idx": mol_idx})
    return result["energy"].detach().cpu().numpy().reshape(-1)


def ef(scores: np.ndarray, labels: np.ndarray, fraction: float = 0.01) -> float:
    n_top = max(1, math.ceil(len(labels) * fraction))
    ranks = rankdata(-scores, method="average")
    return float((labels[ranks <= n_top].sum() / n_top) / labels.mean())


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    ranks = rankdata(scores, method="average")
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-dir", required=True)
    parser.add_argument("--poses", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    args = parser.parse_args()

    target_dir, pose_path, out = Path(args.target_dir), Path(args.poses), Path(args.out)
    pocket_xyz, pocket_z, pocket_charge = read_pocket(target_dir / "pocket.pdb")
    poses = read_poses(pose_path, target_dir / "library.csv.gz")
    print(f"pocket_atoms={len(pocket_z)} poses={len(poses)} cuda={torch.cuda.is_available()}", flush=True)

    started = time.time()
    calc = AIMNet2Calculator("aimnet2-2025", device="cuda")
    print("metadata=" + json.dumps(calc.metadata, default=str, sort_keys=True), flush=True)
    pocket_e = float(eval_systems(calc, [(pocket_xyz, pocket_z, pocket_charge)])[0])
    supported = set((calc.metadata or {}).get("implemented_species") or [])
    results = []
    skipped = []
    for begin in range(0, len(poses), args.batch_size):
        group = poses[begin : begin + args.batch_size]
        valid = []
        systems = []
        for row in group:
            unsupported = sorted(set(map(int, row["numbers"])) - supported) if supported else []
            if unsupported:
                skipped.append({"mol_id": row["mol_id"], "reason": f"unsupported_elements={unsupported}"})
                continue
            complex_xyz = np.concatenate([pocket_xyz, row["coords"]], axis=0)
            complex_z = np.concatenate([pocket_z, row["numbers"]], axis=0)
            ligand_charge = row["charge"]
            complex_charge = pocket_charge + ligand_charge
            systems.extend([(complex_xyz, complex_z, complex_charge), (row["coords"], row["numbers"], ligand_charge)])
            valid.append(row)
        if not valid:
            continue
        energies = eval_systems(calc, systems)
        for i, row in enumerate(valid):
            complex_e, ligand_e = map(float, energies[2 * i : 2 * i + 2])
            interaction_ev = complex_e - pocket_e - ligand_e
            results.append({
                "mol_id": row["mol_id"], "source_index": row["source_index"],
                "label": row["label"], "paff": row["paff"], "formal_charge": row["charge"],
                "smina_affinity_kcal_mol": row["smina_affinity"],
                "aimnet_complex_eV": complex_e, "aimnet_pocket_eV": pocket_e,
                "aimnet_ligand_eV": ligand_e, "aimnet_interaction_eV": interaction_ev,
                "aimnet_interaction_kcal_mol": interaction_ev * EV_TO_KCAL_MOL,
                "aimnet_score": -interaction_ev * EV_TO_KCAL_MOL,
            })
        print(f"scored={len(results)}/{len(poses)}", flush=True)

    frame = pd.DataFrame(results)
    frame["aimnet_rank"] = rankdata(-frame.aimnet_score, method="average")
    frame["smina_rank"] = rankdata(frame.smina_affinity_kcal_mol, method="average")
    frame.sort_values("aimnet_rank").to_csv(out, index=False)
    labels = frame.label.to_numpy(int)
    aimnet_score = frame.aimnet_score.to_numpy(float)
    smina_score = -frame.smina_affinity_kcal_mol.to_numpy(float)
    summary = {
        "model": "aimnet2-2025", "definition": "E_complex - E_pocket - E_ligand",
        "pocket_formal_charge": pocket_charge, "n_input": len(poses), "n_scored": len(frame),
        "n_skipped": len(skipped), "elapsed_seconds": time.time() - started,
        "interaction_kcal_min": float(frame.aimnet_interaction_kcal_mol.min()),
        "interaction_kcal_median": float(frame.aimnet_interaction_kcal_mol.median()),
        "interaction_kcal_max": float(frame.aimnet_interaction_kcal_mol.max()),
        "aimnet_chunk_ef1": ef(aimnet_score, labels), "aimnet_chunk_auc": auc(aimnet_score, labels),
        "smina_chunk_ef1": ef(smina_score, labels), "smina_chunk_auc": auc(smina_score, labels),
        "aimnet_vs_smina_spearman": float(spearmanr(aimnet_score, smina_score).statistic),
        "cpcm_and_strain_terms_included": False,
        "warning": "This active-enriched first chunk is a functionality test, not full-library EF1.",
        "skipped": skipped,
    }
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
