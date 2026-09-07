import argparse
import csv
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

from rdkit import Chem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-csv", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--casf-root", required=True)
    parser.add_argument("--target", default="3ejr")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--output", required=True)
    parser.add_argument("--obabel", required=True)
    args = parser.parse_args()

    out = Path(args.output)
    mol2_dir = out / "mol2"
    sdf_dir = out / "sdf"
    mol2_dir.mkdir(parents=True, exist_ok=True)
    sdf_dir.mkdir(parents=True, exist_ok=True)

    with open(args.top_csv, newline="") as handle:
        rows = list(csv.DictReader(handle))[: args.top_n]
    for i, row in enumerate(rows, 1):
        row["retrieval_chunk_index"] = str(i)

    wanted = {
        f"CASF-2016/decoys_screening/{args.target}/{args.target}_{row['ligand_pdb']}.mol2": row
        for row in rows
    }
    found = set()
    with tarfile.open(args.archive, "r:gz") as tf:
        for member in tf:
            if member.name not in wanted:
                continue
            source = tf.extractfile(member)
            if source is None:
                continue
            destination = mol2_dir / Path(member.name).name
            with source, destination.open("wb") as sink:
                shutil.copyfileobj(source, sink)
            found.add(member.name)
            if len(found) == len(wanted):
                break
    missing = sorted(set(wanted) - found)
    if missing:
        raise FileNotFoundError(missing)

    pocket_candidates = [
        Path(args.casf_root) / "coreset" / args.target / f"{args.target}_pocket.pdb",
        Path(args.casf_root) / "coreset" / args.target / f"{args.target}_protein.pdb",
    ]
    pocket = next((p for p in pocket_candidates if p.exists()), None)
    if pocket is None:
        raise FileNotFoundError(pocket_candidates)
    shutil.copy2(pocket, out / "pocket.pdb")

    combined = out / f"{args.target}_retrieved_top{args.top_n}_allposes.sdf"
    writer = Chem.SDWriter(str(combined))
    pose_counts = {}
    for row in rows:
        ligand = row["ligand_pdb"]
        mol2 = mol2_dir / f"{args.target}_{ligand}.mol2"
        sdf = sdf_dir / f"{args.target}_{ligand}.sdf"
        subprocess.run([args.obabel, "-imol2", str(mol2), "-osdf", "-O", str(sdf)], check=True)
        count = 0
        for pose_index, mol in enumerate(Chem.SDMolSupplier(str(sdf), removeHs=False), 1):
            if mol is None:
                continue
            count += 1
            mol_id = int(row["retrieval_chunk_index"]) * 1000 + pose_index
            mol.SetProp("_Name", f"{ligand}:pose{pose_index}")
            mol.SetIntProp("mol_id", mol_id)
            mol.SetProp("ligand_pdb", ligand)
            mol.SetIntProp("pose_index", pose_index)
            mol.SetIntProp("label", int(row["label"]))
            mol.SetIntProp("ligunity_rank", int(row["ligunity_rank"]))
            mol.SetDoubleProp("ligunity_score", float(row["ligunity_score"]))
            writer.write(mol)
        pose_counts[ligand] = count
    writer.close()
    (out / "preparation.json").write_text(json.dumps({
        "target": args.target,
        "top_n": args.top_n,
        "pocket_source": str(pocket),
        "pose_counts": pose_counts,
        "total_poses": sum(pose_counts.values()),
        "combined_sdf": str(combined),
    }, indent=2) + "\n")
    print(json.dumps({"target": args.target, "pose_counts": pose_counts,
                      "total_poses": sum(pose_counts.values())}, indent=2))


if __name__ == "__main__":
    main()
