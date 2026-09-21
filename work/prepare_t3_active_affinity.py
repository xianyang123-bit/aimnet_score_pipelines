#!/usr/bin/env python3
"""Prepare a balanced active-only T3 affinity-ranking subset for docking."""

from __future__ import annotations

import argparse
import json
import os
import pickle
import subprocess
import shutil
from scipy.spatial import cKDTree
from prepared_pocket import read_pocket
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd


def atom_element(atom_name: str) -> str:
    name = "".join(c for c in str(atom_name).upper() if c.isalpha())
    if name.startswith("CL"):
        return "Cl"
    if name.startswith("BR"):
        return "Br"
    return name[0] if name else "C"


def load_pocket(path: Path) -> dict:
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        value = next(iter(txn.cursor()))[1]
    env.close()
    return pickle.loads(value)


def write_pocket(record: dict, path: Path) -> np.ndarray:
    coords = np.asarray(record["pocket_coordinates"], dtype=float)
    names = record["pocket_atoms"]
    with path.open("w", newline="\n") as handle:
        for i, (name, xyz) in enumerate(zip(names, coords), 1):
            element = atom_element(name)
            handle.write(
                f"ATOM  {i:5d} {str(name)[:4]:<4s} POC A   1    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n"
            )
        handle.write("END\n")
    return coords


def write_json(value: object, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared-pocket-root", required=True, help="Root containing LAYER/UNIPROT/prepared.pdb and .prep.json")
    ap.add_argument("--data", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--env-bin", required=True)
    ap.add_argument("--targets-per-layer", type=int, default=24)
    ap.add_argument("--ligands-per-target", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260904)
    ap.add_argument("--padding", type=float, default=4.0)
    args = ap.parse_args()

    data, work = Path(args.data), Path(args.work)
    if work.exists() and any(work.iterdir()):
        raise FileExistsError("Refusing to reuse a nonempty preparation directory; choose a new --work path")
    work.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(data / "targets.csv.gz")
    actives = pd.read_csv(data / "actives.csv.gz")
    molecules = pd.read_csv(data / "molecules.csv.gz", dtype={"mol_id": np.int64}).set_index("mol_id")
    rng = np.random.default_rng(args.seed)
    task_rows, selected_rows = [], []

    for layer in ("L1", "L2", "L3", "L4"):
        candidates = []
        for row in targets[targets.layer == layer].sort_values("uniprot").itertuples(index=False):
            pocket = data / "pockets_6A" / layer / row.uniprot / f"{row.uniprot}_pocket.lmdb"
            n = int(((actives.layer == layer) & (actives.uniprot == row.uniprot)).sum())
            if pocket.exists() and n >= args.ligands_per_target:
                candidates.append((str(row.uniprot), pocket, n))
        if len(candidates) < args.targets_per_layer:
            raise RuntimeError(f"{layer}: only {len(candidates)} eligible targets")
        chosen_idx = sorted(rng.choice(len(candidates), args.targets_per_layer, replace=False).tolist())
        for idx in chosen_idx:
            up, pocket_lmdb, available = candidates[idx]
            target_dir = work / layer / up
            target_dir.mkdir(parents=True, exist_ok=True)
            subset = actives[(actives.layer == layer) & (actives.uniprot == up)].copy()
            subset = subset.sort_values("mol_id")
            pick = sorted(rng.choice(len(subset), args.ligands_per_target, replace=False).tolist())
            subset = subset.iloc[pick].copy().sort_values("mol_id")
            subset["smiles"] = subset.mol_id.astype(np.int64).map(molecules.smiles)
            if subset.smiles.isna().any():
                raise RuntimeError(f"{layer}/{up}: missing SMILES")
            subset.insert(0, "source_index", np.arange(len(subset), dtype=np.int64))
            subset.to_csv(target_dir / "actives.csv", index=False)

            record = load_pocket(pocket_lmdb)
            coords = write_pocket(record, target_dir / "pocket.reference.pdb")
            prepared = Path(args.prepared_pocket_root) / layer / up / "prepared.pdb"
            prepared_xyz, prepared_z, _ = read_pocket(prepared)
            distances, _ = cKDTree(prepared_xyz[prepared_z != 1]).query(coords)
            if not np.all(distances < 0.02):
                raise ValueError(f"{layer}/{up}: prepared pocket does not retain reference coordinates")
            if (target_dir / "pocket.pdb").exists():
                raise FileExistsError("Use a new work directory for the prepared-pocket protocol")
            shutil.copyfile(prepared, target_dir / "pocket.pdb")
            shutil.copyfile(prepared.with_suffix(".prep.json"), target_dir / "pocket.prep.json")
            lo, hi = coords.min(axis=0) - args.padding, coords.max(axis=0) + args.padding
            center, size = (lo + hi) / 2.0, hi - lo
            box = {"center": center.tolist(), "size": size.tolist(), "padding": args.padding}
            write_json(box, target_dir / "box.json")
            command = [
                str(Path(args.env_bin) / "obabel"), str(target_dir / "pocket.pdb"),
                "-O", str(target_dir / "pocket.pdbqt"), "-xr", "-p", "7.4",
            ]
            proc = subprocess.run(command, text=True, capture_output=True)
            (target_dir / "obabel.log").write_text(proc.stdout + "\n" + proc.stderr)
            if proc.returncode or not (target_dir / "pocket.pdbqt").exists():
                raise RuntimeError(f"{layer}/{up}: receptor conversion failed")
            task_id = len(task_rows)
            task_rows.append({"task_id": task_id, "layer": layer, "uniprot": up,
                              "target_dir": str(target_dir), "n_ligands": len(subset)})
            selected_rows.append({"layer": layer, "uniprot": up,
                                  "available_actives": available, "selected_actives": len(subset),
                                  "pocket_atoms": len(coords), **box})

    pd.DataFrame(task_rows).to_csv(work / "tasks.tsv", sep="\t", index=False)
    pd.DataFrame(selected_rows).to_csv(work / "selected_targets.csv", index=False)
    manifest = {
        "layers": ["L1", "L2", "L3", "L4"],
        "targets_per_layer": args.targets_per_layer,
        "ligands_per_target": args.ligands_per_target,
        "n_targets": len(task_rows),
        "n_complexes": int(sum(x["n_ligands"] for x in task_rows)),
        "seed": args.seed,
        "docking": {"program": "smina", "exhaustiveness": 8, "num_modes": 1,
                    "seed": 1, "box": "6A pocket bounding box + 4A padding"},
        "pocket_protocol": "Validated source-derived, protonated, capped pockets; LMDB is a coordinate reference only.",
    }
    write_json(manifest, work / "manifest.json")
    print(json.dumps(manifest, indent=2))
    print(pd.DataFrame(selected_rows).groupby("layer").agg(
        targets=("uniprot", "size"), complexes=("selected_actives", "sum"),
        mean_pocket_atoms=("pocket_atoms", "mean")).to_string())


if __name__ == "__main__":
    main()

