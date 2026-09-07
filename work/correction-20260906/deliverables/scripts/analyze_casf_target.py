import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def metrics(frame, score, ascending):
    ranked = frame.sort_values(score, ascending=ascending, kind="stable").reset_index(drop=True)
    y = ranked.label.astype(int).to_numpy()
    n = len(y); n_active = int(y.sum()); n1 = max(1, math.ceil(0.01 * n))
    return {"n": n, "n_active": n_active, "ef1_top_count": n1,
            "ef1": float((y[:n1].sum() / n1) / (n_active / n)) if n_active else None,
            "active_top1": int(y[:1].sum()), "active_top5": int(y[:5].sum()),
            "active_top10": int(y[:10].sum()),
            "active_ranks": {row.ligand_pdb: i + 1 for i, row in ranked.iterrows() if int(row.label) == 1}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-dir", required=True)
    parser.add_argument("--target", required=True)
    args = parser.parse_args()
    root = Path(args.target_dir)
    retrieval = pd.read_csv(root / "ligunity_all285.csv")
    top20 = retrieval.nsmallest(20, "ligunity_rank")
    interaction = pd.read_csv(root / "aimnet_interaction_bestpose.csv")
    composite = pd.read_csv(root / "aimnet2_composite_rerank.csv")
    composite = composite[composite.aimnet2_composite_kcal_mol.notna()].copy()
    summary = {
        "target": args.target,
        "full_285_ligunity": metrics(retrieval, "ligunity_rank", True),
        "top20_ligunity": metrics(top20, "ligunity_rank", True),
        "top20_aimnet_interaction": metrics(interaction, "aimnet_interaction_kcal_mol", True),
        "top20_aimnet2_composite": metrics(composite, "aimnet2_composite_kcal_mol", True),
    }
    merged = top20[["ligand_pdb", "smiles", "label", "ligunity_score", "ligunity_rank"]].merge(
        interaction[["ligand_pdb", "pose_index", "aimnet_interaction_kcal_mol", "aimnet_interaction_rank"]],
        on="ligand_pdb", how="left").merge(
        composite[["ligand_pdb", "desolvation_kcal_mol", "local_strain_kcal_mol",
                   "aimnet2_composite_kcal_mol", "aimnet2_score_rank", "optimization_converged"]],
        on="ligand_pdb", how="left")
    summary["spearman_ligunity_vs_interaction"] = float(spearmanr(
        merged.ligunity_rank, merged.aimnet_interaction_kcal_mol).statistic)
    summary["spearman_ligunity_vs_composite"] = float(spearmanr(
        merged.ligunity_rank, merged.aimnet2_composite_kcal_mol).statistic)
    summary["n_composite_converged"] = int(merged.optimization_converged.fillna(False).sum())
    merged.sort_values("aimnet2_score_rank").to_csv(root / "ranking_comparison_all20.csv", index=False)
    (root / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
