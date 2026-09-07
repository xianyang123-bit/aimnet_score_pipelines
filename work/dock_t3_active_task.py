#!/usr/bin/env python3
"""Dock one balanced T3 active-only target task with smina."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


def write_json(value: object, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--env-bin", required=True)
    ap.add_argument("--task-id", type=int)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--timeout-seconds", type=int, default=600)
    args = ap.parse_args()
    task_id = args.task_id if args.task_id is not None else int(os.environ["SLURM_ARRAY_TASK_ID"])
    work = Path(args.work)
    task = pd.read_csv(work / "tasks.tsv", sep="\t").iloc[task_id]
    target_dir = Path(task.target_dir)
    output = target_dir / "poses.sdf"
    status_path = target_dir / "dock_status.json"
    if status_path.exists() and output.exists():
        old = json.loads(status_path.read_text())
        if old.get("complete") and old.get("atom_count_mismatches", 1) == 0:
            print(f"SKIP {task.layer}/{task.uniprot}")
            return

    frame = pd.read_csv(target_dir / "actives.csv")
    input_sdf = target_dir / "input_3d.sdf"
    writer = Chem.SDWriter(str(input_sdf))
    failures = []
    input_molecules = []
    input_atom_counts = {}
    for row in frame.itertuples(index=False):
        mol = Chem.MolFromSmiles(str(row.smiles))
        if mol is None:
            failures.append({"mol_id": int(row.mol_id), "stage": "parse"})
            continue
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(int(row.mol_id) % 2147483646) + 1
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                failures.append({"mol_id": int(row.mol_id), "stage": "embed"})
                continue
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
        mol.SetProp("_Name", f"{task.layer}:{task.uniprot}:{int(row.mol_id)}")
        mol.SetIntProp("mol_id", int(row.mol_id))
        mol.SetIntProp("source_index", int(row.source_index))
        mol.SetProp("layer", str(task.layer))
        mol.SetProp("uniprot", str(task.uniprot))
        mol.SetDoubleProp("paff", float(row.paff))
        mol.SetProp("smiles", str(row.smiles))
        input_atom_counts[int(row.mol_id)] = mol.GetNumAtoms()
        input_molecules.append(Chem.Mol(mol))
        writer.write(mol)
    writer.close()

    box = json.loads((target_dir / "box.json").read_text())
    command = [
        str(Path(args.env_bin) / "smina"), "-r", str(target_dir / "pocket.pdbqt"),
        "-l", str(input_sdf), "--center_x", str(box["center"][0]),
        "--center_y", str(box["center"][1]), "--center_z", str(box["center"][2]),
        "--size_x", str(box["size"][0]), "--size_y", str(box["size"][1]),
        "--size_z", str(box["size"][2]), "--exhaustiveness", str(args.exhaustiveness),
        "--num_modes", "1", "--cpu", str(args.cpus), "--seed", "1", "-o", str(output),
    ]
    # GNU timeout terminates the whole docking command even when smina leaves
    # worker processes holding stdout/stderr open.  The outer Python timeout is
    # only a final guard around timeout's TERM/KILL sequence.
    guarded_command = [
        "timeout", "--kill-after=15s", f"{args.timeout_seconds}s", *command
    ]
    try:
        proc = subprocess.run(
            guarded_command,
            text=True,
            capture_output=True,
            timeout=args.timeout_seconds + 45,
        )
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        proc = subprocess.CompletedProcess(guarded_command, 124, stdout, stderr + "\nTIMEOUT\n")
        timed_out = True
    if proc.returncode in (124, 137):
        timed_out = True
    (target_dir / "smina.log").write_text("COMMAND: " + " ".join(guarded_command) + "\n\n" + proc.stdout + "\n" + proc.stderr)
    raw_poses = ([m for m in Chem.SDMolSupplier(
        str(output), removeHs=False, sanitize=False, strictParsing=False
    ) if m] if output.exists() else [])
    hydrogenated = target_dir / "poses.hydrogenated.sdf"
    h_command = [str(Path(args.env_bin) / "obabel"), str(output), "-O", str(hydrogenated), "-h"]
    h_proc = subprocess.run(h_command, text=True, capture_output=True)
    (target_dir / "obabel_ligand_h.log").write_text(h_proc.stdout + "\n" + h_proc.stderr)
    poses = [m for m in Chem.SDMolSupplier(str(hydrogenated), removeHs=False) if m] if hydrogenated.exists() else []
    reconstruction_failures = []
    # Smina preserves geometry and its affinity field but drops arbitrary SDF
    # properties. Restore identity/labels from the one-pose-per-input ordering;
    # prefer the title when it survived unchanged.
    input_by_name = {m.GetProp("_Name"): m for m in input_molecules}
    for i, mol in enumerate(poses):
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        src = input_by_name.get(name, input_molecules[i] if i < len(input_molecules) else None)
        if src is not None:
            for key in ("_Name", "mol_id", "source_index", "layer", "uniprot", "paff", "smiles"):
                if src.HasProp(key):
                    mol.SetProp(key, src.GetProp(key))
    repaired = target_dir / "poses.repaired.sdf"
    repaired_writer = Chem.SDWriter(str(repaired))
    for mol in poses:
        repaired_writer.write(mol)
    repaired_writer.close()
    os.replace(repaired, output)

    missing_props, atom_mismatches = 0, 0
    for mol in poses:
        if not mol.HasProp("mol_id"):
            missing_props += 1
            continue
        mid = int(mol.GetProp("mol_id"))
        if mid in input_atom_counts and mol.GetNumAtoms() != input_atom_counts[mid]:
            atom_mismatches += 1
    status = {
        "task_id": int(task_id), "layer": str(task.layer), "uniprot": str(task.uniprot),
        "requested": len(frame), "embedded": len(input_atom_counts), "output_poses": len(poses),
        "embed_failures": len(failures), "smina_returncode": proc.returncode,
        "timed_out": timed_out,
        "missing_mol_id_properties": missing_props, "atom_count_mismatches": atom_mismatches,
        "raw_output_poses": len(raw_poses), "reconstruction_failures": len(reconstruction_failures),
        "hydrogenation_returncode": h_proc.returncode,
        "complete": proc.returncode == 0 and len(poses) == len(input_atom_counts) and missing_props == 0 and atom_mismatches == 0,
    }
    all_failures = failures + reconstruction_failures
    if all_failures:
        pd.DataFrame(all_failures).to_csv(target_dir / "dock_failures.csv", index=False)
    write_json(status, status_path)
    print(json.dumps(status, indent=2))
    if not status["complete"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
