import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/data/user_data/xianyang/casf-2016/ligunity")
RERANK = ROOT / "rerank-3ebp-top20"


def metrics(frame, score, ascending, label="label"):
    ranked = frame.sort_values(score, ascending=ascending, kind="stable").reset_index(drop=True)
    y = ranked[label].astype(int).to_numpy()
    n = len(y); n_active = int(y.sum())
    n1 = max(1, math.ceil(0.01 * n))
    ef1 = float((y[:n1].sum() / n1) / (n_active / n)) if n_active else None
    return {
        "n": n, "n_active": n_active, "ef1_top_count": n1, "ef1": ef1,
        "active_top1": int(y[:1].sum()), "active_top5": int(y[:5].sum()),
        "active_top10": int(y[:10].sum()),
        "active_ranks": {row.ligand_pdb: i + 1 for i, row in ranked.iterrows() if int(row[label]) == 1},
    }


def main():
    retrieval = pd.read_csv(ROOT / "results/3ebp_ligunity_all285.csv")
    top20 = retrieval.nsmallest(20, "ligunity_rank")
    interaction = pd.read_csv(RERANK / "aimnet_interaction_bestpose.csv")
    composite = pd.read_csv(RERANK / "aimnet2_composite_rerank.csv")
    composite = composite[composite.aimnet2_composite_kcal_mol.notna()].copy()
    summary = {
        "target": "3ebp",
        "target_selection": "first official CASF screening target whose pocket has only public AIMNet2-supported elements",
        "full_285_ligunity": metrics(retrieval, "ligunity_rank", True),
        "top20_ligunity": metrics(top20, "ligunity_rank", True),
        "top20_aimnet_interaction_best_of_100_poses": metrics(interaction, "aimnet_interaction_kcal_mol", True),
        "top20_aimnet2_composite": metrics(composite, "aimnet2_composite_kcal_mol", True),
        "composite_reference_warning": "Public-checkpoint reconstruction of Eint + Edesolv + ELCSE; not the unreleased trained AIMNet2(Score) checkpoint.",
        "pose_protocol": "100 provided CASF docked poses per retrieved ligand; select minimum AIMNet2 interaction energy, then compute composite score on that pose.",
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
    merged.sort_values("aimnet2_score_rank").to_csv(RERANK / "ranking_comparison_all20.csv", index=False)
    (RERANK / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(merged.sort_values("aimnet2_score_rank").to_string(index=False))


if __name__ == "__main__":
    main()
