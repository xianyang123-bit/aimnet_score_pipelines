import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

SUPPORTED = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 53}


def pocket_elements(path: Path):
    table = Chem.GetPeriodicTable()
    numbers = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        symbol = line[76:78].strip()
        if not symbol:
            name = "".join(c for c in line[12:16] if c.isalpha())
            symbol = name[0] if name else "C"
        numbers.append(table.GetAtomicNumber(symbol.title()))
    return numbers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--casf-root", required=True)
    parser.add_argument("--results", required=True)
    args = parser.parse_args()
    root = Path(args.casf_root)
    results = Path(args.results)
    ids = json.loads((results / "PDBBind/test_pdbbind_ids.json").read_text())
    mol = np.load(results / "PDBBind/test_mol_reps.npy")
    prot = np.load(results / "PDBBind/test_pocket_reps.npy")
    smis = json.loads((results / "PDBBind/test_mol_smis.json").read_text())
    rows = []
    target_lines = (root / "power_screening/TargetInfo.dat").read_text().splitlines()[1:]
    groups = [line.split() for line in target_lines if line.strip()]
    for group in groups:
        target, actives = group[0], group[1:]
        pocket = root / "coreset" / target / f"{target}_pocket.pdb"
        if not pocket.exists():
            pocket = root / "coreset" / target / f"{target}_protein.pdb"
        z = pocket_elements(pocket)
        unsupported = sorted(set(z) - SUPPORTED)
        score = prot[ids.index(target)] @ mol.T
        labels = np.array([int(x in set(actives)) for x in ids])
        order = np.argsort(-score, kind="stable")
        rank = np.empty(len(order), dtype=int); rank[order] = np.arange(1, len(order) + 1)
        n_top = max(1, math.ceil(0.01 * len(ids)))
        ef1 = (labels[order[:n_top]].sum() / n_top) / labels.mean()
        rows.append({
            "target": target, "pocket": str(pocket), "pocket_atoms": len(z),
            "unsupported_atomic_numbers": ",".join(map(str, unsupported)),
            "n_actives": int(labels.sum()), "actives": ",".join(actives),
            "ef1": float(ef1), "actives_top20": int(labels[order[:20]].sum()),
            "active_ranks": ",".join(f"{x}:{rank[ids.index(x)]}" for x in actives),
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(results / "casf_target_support_and_retrieval.csv", index=False)
    selected = frame[frame.unsupported_atomic_numbers == ""].iloc[0]
    target = selected.target
    actives = selected.actives.split(",")
    score = prot[ids.index(target)] @ mol.T
    labels = np.array([int(x in set(actives)) for x in ids])
    order = np.argsort(-score, kind="stable")
    rank = np.empty(len(order), dtype=int); rank[order] = np.arange(1, len(order) + 1)
    ranking = pd.DataFrame({"ligand_pdb": ids, "smiles": smis, "ligunity_score": score,
                            "ligunity_rank": rank, "label": labels}).sort_values("ligunity_rank")
    ranking.to_csv(results / f"{target}_ligunity_all285.csv", index=False)
    ranking.head(50).to_csv(results / f"{target}_ligunity_top50.csv", index=False)
    summary = {"selection_rule": "first CASF TargetInfo target with AIMNet2-supported pocket elements",
               **selected.to_dict()}
    (results / "selected_supported_target.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(frame.head(15).to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(ranking.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
