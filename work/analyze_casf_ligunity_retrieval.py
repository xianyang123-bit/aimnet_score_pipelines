import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--target", default="3ejr")
    parser.add_argument("--actives", nargs="+", required=True)
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    results = Path(args.results)
    mol = np.load(results / "PDBBind/test_mol_reps.npy")
    prot = np.load(results / "PDBBind/test_pocket_reps.npy")
    ids = json.loads((results / "PDBBind/test_pdbbind_ids.json").read_text())
    smis = json.loads((results / "PDBBind/test_mol_smis.json").read_text())
    if len(mol) != len(ids) or len(prot) != len(ids):
        raise ValueError((mol.shape, prot.shape, len(ids)))
    target_index = ids.index(args.target)
    scores = prot[target_index] @ mol.T
    active_set = set(args.actives)
    labels = np.array([int(pdb in active_set) for pdb in ids], dtype=int)
    order = np.argsort(-scores, kind="stable")
    rank = np.empty(len(order), dtype=int)
    rank[order] = np.arange(1, len(order) + 1)
    top_n = min(args.top_n, len(ids))
    frame = pd.DataFrame({
        "ligand_pdb": ids,
        "smiles": smis,
        "ligunity_score": scores,
        "ligunity_rank": rank,
        "label": labels,
    }).sort_values("ligunity_rank")
    frame.to_csv(results / f"{args.target}_ligunity_all285.csv", index=False)
    frame.head(top_n).to_csv(results / f"{args.target}_ligunity_top{top_n}.csv", index=False)

    n_top = max(1, math.ceil(0.01 * len(labels)))
    active_top = int(labels[order[:n_top]].sum())
    ef1 = (active_top / n_top) / labels.mean()
    active_ranks = {pdb: int(rank[i]) for i, pdb in enumerate(ids) if labels[i]}
    summary = {
        "dataset": "CASF-2016 screening library",
        "retrieval_model": "LigUnity protein_ranking_vs checkpoint_avg_41-50.pt",
        "target": args.target,
        "n_candidates": len(ids),
        "n_actives": int(labels.sum()),
        "top_n_exported": top_n,
        "actives_in_top_n": int(labels[order[:top_n]].sum()),
        "ef1_top_count": n_top,
        "ef1_actives_in_top": active_top,
        "ef1": float(ef1),
        "active_ranks": active_ranks,
    }
    (results / f"{args.target}_ligunity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(frame.head(top_n).to_string(index=False))


if __name__ == "__main__":
    main()
