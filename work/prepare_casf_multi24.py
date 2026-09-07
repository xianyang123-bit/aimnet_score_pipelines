import argparse
import csv
import json
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--casf-root", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--retrieval-results", required=True)
    parser.add_argument("--support-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--obabel", required=True)
    parser.add_argument("--n-targets", type=int, default=24)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    casf = Path(args.casf_root)
    retrieval = Path(args.retrieval_results)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ids = json.loads((retrieval / "PDBBind/test_pdbbind_ids.json").read_text())
    smis = json.loads((retrieval / "PDBBind/test_mol_smis.json").read_text())
    mol = np.load(retrieval / "PDBBind/test_mol_reps.npy")
    prot = np.load(retrieval / "PDBBind/test_pocket_reps.npy")
    support = pd.read_csv(args.support_csv).fillna("")
    selected = support[support.unsupported_atomic_numbers.astype(str) == ""].head(args.n_targets).copy()
    if len(selected) != args.n_targets:
        raise ValueError(f"only {len(selected)} supported targets")
    selected.insert(0, "array_index", np.arange(len(selected)))
    selected.to_csv(out / "selected_targets.csv", index=False)

    wanted = {}
    target_rows = {}
    for selected_row in selected.itertuples(index=False):
        target = selected_row.target
        active_set = set(str(selected_row.actives).split(","))
        score = prot[ids.index(target)] @ mol.T
        order = np.argsort(-score, kind="stable")
        rank = np.empty(len(order), dtype=int)
        rank[order] = np.arange(1, len(order) + 1)
        frame = pd.DataFrame({"ligand_pdb": ids, "smiles": smis, "ligunity_score": score,
                              "ligunity_rank": rank,
                              "label": [int(x in active_set) for x in ids]}).sort_values("ligunity_rank")
        target_dir = out / target
        target_dir.mkdir(exist_ok=True)
        frame.to_csv(target_dir / "ligunity_all285.csv", index=False)
        top = frame.head(args.top_n).copy()
        top.to_csv(target_dir / f"ligunity_top{args.top_n}.csv", index=False)
        target_rows[target] = top.to_dict("records")
        for row in target_rows[target]:
            name = f"CASF-2016/decoys_screening/{target}/{target}_{row['ligand_pdb']}.mol2"
            wanted[name] = (target, row)

    tmp_root = Path(tempfile.mkdtemp(prefix="casf-multi24-", dir=out))
    found = set()
    try:
        with tarfile.open(args.archive, "r:gz") as tf:
            for member in tf:
                item = wanted.get(member.name)
                if item is None:
                    continue
                target, row = item
                source = tf.extractfile(member)
                if source is None:
                    continue
                destination = tmp_root / target / Path(member.name).name
                destination.parent.mkdir(parents=True, exist_ok=True)
                with source, destination.open("wb") as sink:
                    shutil.copyfileobj(source, sink)
                found.add(member.name)
                if len(found) % 50 == 0:
                    print(f"extracted={len(found)}/{len(wanted)}", flush=True)
                if len(found) == len(wanted):
                    break
        missing = sorted(set(wanted) - found)
        if missing:
            raise FileNotFoundError(missing[:20])

        for target_index, selected_row in enumerate(selected.itertuples(index=False), 1):
            target = selected_row.target
            target_dir = out / target
            pocket = casf / "coreset" / target / f"{target}_pocket.pdb"
            if not pocket.exists():
                pocket = casf / "coreset" / target / f"{target}_protein.pdb"
            shutil.copy2(pocket, target_dir / "pocket.pdb")
            combined = target_dir / f"{target}_top{args.top_n}_allposes.sdf"
            writer = Chem.SDWriter(str(combined))
            counts = {}
            for chunk_index, row in enumerate(target_rows[target], 1):
                ligand = row["ligand_pdb"]
                mol2 = tmp_root / target / f"{target}_{ligand}.mol2"
                sdf = tmp_root / target / f"{target}_{ligand}.sdf"
                completed = subprocess.run(
                    [args.obabel, "-imol2", str(mol2), "-osdf", "-O", str(sdf)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                )
                if completed.returncode != 0:
                    raise RuntimeError(completed.stdout[-2000:])
                count = 0
                for pose_index, pose in enumerate(Chem.SDMolSupplier(str(sdf), removeHs=False), 1):
                    if pose is None:
                        continue
                    count += 1
                    pose.SetProp("_Name", f"{ligand}:pose{pose_index}")
                    pose.SetIntProp("mol_id", chunk_index * 1000 + pose_index)
                    pose.SetProp("ligand_pdb", ligand)
                    pose.SetIntProp("pose_index", pose_index)
                    pose.SetIntProp("label", int(row["label"]))
                    pose.SetIntProp("ligunity_rank", int(row["ligunity_rank"]))
                    pose.SetDoubleProp("ligunity_score", float(row["ligunity_score"]))
                    writer.write(pose)
                counts[ligand] = count
            writer.close()
            prep = {"target": target, "top_n": args.top_n, "pocket_source": str(pocket),
                    "pose_counts": counts, "total_poses": sum(counts.values()),
                    "combined_sdf": str(combined)}
            (target_dir / "preparation.json").write_text(json.dumps(prep, indent=2) + "\n")
            print(f"prepared target={target} ({target_index}/{len(selected)}) poses={sum(counts.values())}", flush=True)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
