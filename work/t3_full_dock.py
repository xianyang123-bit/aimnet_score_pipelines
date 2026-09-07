#!/usr/bin/env python3
"""Resumable full-library smina docking for the T3 L4 benchmark subset.

The target set and docking settings reproduce the existing T6 smina experiment,
but every active and decoy is docked instead of only a retrieval top-200.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import pickle
import shutil
import subprocess
import tempfile
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from scipy.stats import rankdata, spearmanr

RDLogger.DisableLog("rdApp.*")

DEFAULT_DATA = Path("/data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset")
DEFAULT_WORK = Path("/data/user_data/xianyang/t3-aimnet-full/full-dock-L4")
DEFAULT_ENV = Path("/home/xianyang/miniconda3/envs/t3-aimnet/bin")

# Same 20 high-quality L4 targets selected by physics/prep_dock.py for T6.
TARGETS = [
    "B4RNM9", "I6WXK4", "J9SQF3", "O00187", "O00411",
    "O00443", "O00748", "O00750", "O00767", "O14578",
    "O14874", "O14929", "O14980", "O15393", "O42275",
    "O54890", "O60427", "O75106", "O75143", "O75365",
]


def atomic_element(atom_name: str) -> str:
    """Infer a protein heavy-atom element from its PDB atom name."""
    name = "".join(c for c in str(atom_name).upper() if c.isalpha())
    if not name:
        return "C"
    # The shipped pockets are standard protein atoms. CA is alpha carbon, not calcium.
    if name.startswith("CL"):
        return "Cl"
    if name.startswith("BR"):
        return "Br"
    return name[0]


def load_pocket(path: Path) -> dict:
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        value = next(iter(txn.cursor()))[1]
    env.close()
    return pickle.loads(value)


def write_pocket_pdb(record: dict, path: Path) -> np.ndarray:
    coords = np.asarray(record["pocket_coordinates"], dtype=float)
    atom_names = record["pocket_atoms"]
    if coords.shape != (len(atom_names), 3):
        raise ValueError(f"bad pocket shape: {coords.shape} for {len(atom_names)} atoms")
    with path.open("w", newline="\n") as handle:
        for i, (name, xyz) in enumerate(zip(atom_names, coords), 1):
            element = atomic_element(name)
            handle.write(
                f"ATOM  {i:5d} {str(name)[:4]:<4s} POC A   1    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n"
            )
        handle.write("END\n")
    return coords


def dump_json(value: object, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def prepare(args: argparse.Namespace) -> None:
    data, work = Path(args.data), Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(data / "targets.csv.gz")
    got = targets.set_index("uniprot")
    for up in TARGETS:
        if up not in got.index:
            raise ValueError(f"target absent from normalized dataset: {up}")
        if got.loc[up, "layer"] != "L4":
            raise ValueError(f"target {up} is not L4")

    molecules = pd.read_csv(data / "molecules.csv.gz", dtype={"mol_id": np.int64})
    molecules = molecules.set_index("mol_id")
    active = pd.read_csv(data / "actives.csv.gz")
    decoy = pd.read_csv(data / "decoys.csv.gz")
    active = active[(active.layer == "L4") & active.uniprot.isin(TARGETS)]
    decoy = decoy[(decoy.layer == "L4") & decoy.uniprot.isin(TARGETS)]

    task_rows: list[dict] = []
    target_meta: dict[str, dict] = {}
    task_id = 0
    for up in TARGETS:
        target_dir = work / up
        target_dir.mkdir(exist_ok=True)
        (target_dir / "chunks").mkdir(exist_ok=True)
        (target_dir / "logs").mkdir(exist_ok=True)
        (target_dir / "status").mkdir(exist_ok=True)

        a = active[active.uniprot == up][["mol_id", "paff"]].copy()
        a["label"] = 1
        d = decoy[decoy.uniprot == up][["mol_id"]].copy()
        d["paff"] = np.nan
        d["label"] = 0
        library = pd.concat([a, d], ignore_index=True)
        library["mol_id"] = library.mol_id.astype(np.int64)
        library["smiles"] = library.mol_id.map(molecules.smiles)
        if library.smiles.isna().any():
            raise ValueError(f"{up}: missing {int(library.smiles.isna().sum())} SMILES")
        if library.mol_id.duplicated().any():
            raise ValueError(f"{up}: duplicate molecule pairs")
        library.insert(0, "source_index", np.arange(len(library), dtype=np.int64))
        library.to_csv(target_dir / "library.csv.gz", index=False)

        pocket_path = data / "pockets_6A" / "L4" / up / f"{up}_pocket.lmdb"
        if not pocket_path.exists():
            raise FileNotFoundError(pocket_path)
        record = load_pocket(pocket_path)
        coords = write_pocket_pdb(record, target_dir / "pocket.pdb")
        lo, hi = coords.min(axis=0) - args.padding, coords.max(axis=0) + args.padding
        center, size = (lo + hi) / 2.0, hi - lo
        command = [
            str(Path(args.env_bin) / "obabel"), str(target_dir / "pocket.pdb"),
            "-O", str(target_dir / "pocket.pdbqt"), "-xr", "-p", "7.4",
        ]
        proc = subprocess.run(command, text=True, capture_output=True)
        (target_dir / "obabel.log").write_text(proc.stdout + "\n" + proc.stderr)
        if proc.returncode or not (target_dir / "pocket.pdbqt").exists():
            raise RuntimeError(f"{up}: receptor conversion failed: {proc.stderr[-500:]}")

        box = {"center": center.tolist(), "size": size.tolist(), "padding": args.padding}
        dump_json(box, target_dir / "box.json")
        target_meta[up] = {
            "n": int(len(library)), "n_actives": int(library.label.sum()),
            "n_decoys": int((library.label == 0).sum()), "pocket_atoms": int(len(coords)),
            **box,
        }
        for start in range(0, len(library), args.chunk_size):
            stop = min(start + args.chunk_size, len(library))
            task_rows.append({
                "task_id": task_id, "uniprot": up, "start": start,
                "stop": stop, "count": stop - start,
            })
            task_id += 1

    with (work / "tasks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=task_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(task_rows)
    dump_json({
        "layer": "L4", "targets": TARGETS, "chunk_size": args.chunk_size,
        "n_tasks": len(task_rows), "n_ligands": int(sum(x["count"] for x in task_rows)),
        "docking": {"program": "smina", "exhaustiveness": 8, "num_modes": 1, "seed": 1},
        "target_metadata": target_meta,
    }, work / "manifest.json")
    print(json.dumps(json.loads((work / "manifest.json").read_text()), indent=2))


def embed_ligand(row: pd.Series) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(str(row.smiles))
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(int(row.mol_id) % 2147483646) + 1
    params.useRandomCoords = False
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    mol.SetProp("_Name", f"{row.uniprot}:{int(row.mol_id)}")
    mol.SetIntProp("mol_id", int(row.mol_id))
    mol.SetIntProp("source_index", int(row.source_index))
    mol.SetIntProp("label", int(row.label))
    if not pd.isna(row.paff):
        mol.SetDoubleProp("paff", float(row.paff))
    return mol


def count_sdf(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(m is not None for m in Chem.SDMolSupplier(str(path), removeHs=False))


def dock_task(args: argparse.Namespace) -> None:
    work = Path(args.work)
    task_id = args.task_id
    if task_id is None:
        task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    tasks = pd.read_csv(work / "tasks.tsv", sep="\t")
    if task_id < 0 or task_id >= len(tasks):
        raise IndexError(f"task {task_id} outside 0..{len(tasks)-1}")
    task = tasks.iloc[task_id]
    up, start, stop = str(task.uniprot), int(task.start), int(task.stop)
    target_dir = work / up
    stem = f"{start:06d}_{stop:06d}"
    status_path = target_dir / "status" / f"{stem}.json"
    pose_path = target_dir / "chunks" / f"poses_{stem}.sdf.gz"
    if status_path.exists():
        old = json.loads(status_path.read_text())
        if old.get("complete") and pose_path.exists():
            print(f"SKIP complete {up} {stem}")
            return

    library = pd.read_csv(target_dir / "library.csv.gz")
    frame = library.iloc[start:stop].copy()
    frame["uniprot"] = up
    scratch_root = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR")
    tmp_dir = Path(tempfile.mkdtemp(prefix=f"t3dock-{up}-{stem}-", dir=scratch_root))
    input_sdf, output_sdf = tmp_dir / "input.sdf", tmp_dir / "poses.sdf"
    failures: list[dict] = []
    writer = Chem.SDWriter(str(input_sdf))
    embedded = 0
    for row in frame.itertuples(index=False):
        mol = embed_ligand(pd.Series(row._asdict()))
        if mol is None:
            failures.append({"source_index": int(row.source_index), "mol_id": int(row.mol_id), "reason": "embed"})
            continue
        writer.write(mol)
        embedded += 1
    writer.close()
    box = json.loads((target_dir / "box.json").read_text())
    command = [
        str(Path(args.env_bin) / "smina"), "-r", str(target_dir / "pocket.pdbqt"),
        "-l", str(input_sdf), "--center_x", str(box["center"][0]),
        "--center_y", str(box["center"][1]), "--center_z", str(box["center"][2]),
        "--size_x", str(box["size"][0]), "--size_y", str(box["size"][1]),
        "--size_z", str(box["size"][2]), "--exhaustiveness", str(args.exhaustiveness),
        "--num_modes", "1", "--cpu", str(args.cpus), "--seed", "1", "-o", str(output_sdf),
    ]
    proc = subprocess.run(command, text=True, capture_output=True)
    log_path = target_dir / "logs" / f"smina_{stem}.log"
    log_path.write_text("COMMAND: " + " ".join(command) + "\n\n" + proc.stdout + "\n" + proc.stderr)
    output_count = count_sdf(output_sdf)
    if output_count:
        tmp_gz = pose_path.with_suffix(pose_path.suffix + ".tmp")
        with output_sdf.open("rb") as src, gzip.open(tmp_gz, "wb", compresslevel=3) as dst:
            shutil.copyfileobj(src, dst)
        os.replace(tmp_gz, pose_path)
    if failures:
        pd.DataFrame(failures).to_csv(target_dir / "status" / f"failures_{stem}.csv", index=False)
    status = {
        "task_id": task_id, "uniprot": up, "start": start, "stop": stop,
        "attempted": len(frame), "embedded": embedded, "embed_failures": len(failures),
        "smina_returncode": proc.returncode, "output_count": output_count,
        "complete": proc.returncode == 0 and output_count == embedded,
        "pose_file": str(pose_path),
    }
    dump_json(status, status_path)
    print(json.dumps(status, indent=2))
    shutil.rmtree(tmp_dir, ignore_errors=True)
    if not status["complete"]:
        raise SystemExit(2)


def affinity_from_mol(mol: Chem.Mol) -> float | None:
    props = mol.GetPropsAsDict()
    preferred = ["minimizedAffinity", "affinity", "vina_affinity"]
    for key in preferred:
        if key in props:
            try:
                return float(props[key])
            except (TypeError, ValueError):
                pass
    for key, value in props.items():
        if "affinity" in key.lower():
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    return None


def enrichment_factor(scores: np.ndarray, labels: np.ndarray, fraction: float) -> float:
    n_top = max(1, math.ceil(len(labels) * fraction))
    ranks = rankdata(-scores, method="average")
    hits = int(labels[ranks <= n_top].sum())
    return (hits / n_top) / (labels.sum() / len(labels))


def collect(args: argparse.Namespace) -> None:
    work = Path(args.work)
    manifest = json.loads((work / "manifest.json").read_text())
    summaries = []
    for up in manifest["targets"]:
        target_dir = work / up
        library = pd.read_csv(target_dir / "library.csv.gz")
        scores: dict[int, float] = {}
        prop_names: set[str] = set()
        for path in sorted((target_dir / "chunks").glob("poses_*.sdf.gz")):
            with gzip.open(path, "rb") as handle:
                supplier = Chem.ForwardSDMolSupplier(handle, removeHs=False)
                for mol in supplier:
                    if mol is None:
                        continue
                    prop_names.update(mol.GetPropNames())
                    mol_id = int(mol.GetProp("mol_id")) if mol.HasProp("mol_id") else int(mol.GetProp("_Name").split(":")[-1])
                    affinity = affinity_from_mol(mol)
                    if affinity is not None:
                        scores[mol_id] = affinity
        library["smina_affinity"] = library.mol_id.map(scores)
        library["smina_score"] = -library.smina_affinity
        finite = library.smina_score.notna()
        ranking_score = library.smina_score.fillna(-np.inf).to_numpy(float)
        labels = library.label.to_numpy(int)
        library["smina_rank"] = rankdata(-ranking_score, method="average")
        library.sort_values(["smina_rank", "source_index"]).to_csv(target_dir / "smina_ranking.csv.gz", index=False)
        active_scored = library[(library.label == 1) & finite & library.paff.notna()]
        rho = float("nan")
        if len(active_scored) >= 2:
            rho = float(spearmanr(active_scored.smina_score, active_scored.paff).statistic)
        summaries.append({
            "uniprot": up, "n": len(library), "n_actives": int(labels.sum()),
            "n_scored": int(finite.sum()), "coverage": float(finite.mean()),
            "ef1": enrichment_factor(ranking_score, labels, 0.01),
            "active_affinity_spearman": rho,
            "sdf_properties": ";".join(sorted(prop_names)),
        })
    out = pd.DataFrame(summaries)
    out.to_csv(work / "smina_full_L4_summary.csv", index=False)
    print(out.to_string(index=False))
    print("\nMacro means:\n" + out[["coverage", "ef1", "active_affinity_spearman"]].mean().to_string())


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--data", default=str(DEFAULT_DATA))
    p.add_argument("--work", default=str(DEFAULT_WORK))
    p.add_argument("--env-bin", default=str(DEFAULT_ENV))
    sub = p.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--chunk-size", type=int, default=250)
    prep.add_argument("--padding", type=float, default=4.0)
    prep.set_defaults(func=prepare)
    dock = sub.add_parser("dock-task")
    dock.add_argument("--task-id", type=int)
    dock.add_argument("--cpus", type=int, default=8)
    dock.add_argument("--exhaustiveness", type=int, default=8)
    dock.set_defaults(func=dock_task)
    coll = sub.add_parser("collect")
    coll.set_defaults(func=collect)
    return p


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
