import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/data/user_data/xianyang/casf-2016/ligunity/multi24")


def main():
    targets = pd.read_csv(ROOT / "selected_targets.csv")
    rows = []
    for record in targets.itertuples(index=False):
        target = record.target
        summary = json.loads((ROOT / target / "comparison_summary.json").read_text())
        interaction_run = json.loads((ROOT / target / "aimnet_interaction_bestpose.summary.json").read_text())
        r = summary["top20_ligunity"]
        i = summary["top20_aimnet_interaction"]
        c = summary["top20_aimnet2_composite"]
        rows.append({
            "target": target, "pocket_atoms": record.pocket_atoms,
            "n_input_poses": interaction_run["n_input_poses"],
            "n_scored_poses": interaction_run["n_scored_poses"],
            "n_actives_full": summary["full_285_ligunity"]["n_active"],
            "full_retrieval_ef1": summary["full_285_ligunity"]["ef1"],
            "retrieval_actives_top5": r["active_top5"], "retrieval_actives_top10": r["active_top10"],
            "interaction_ef1": i["ef1"], "interaction_actives_top5": i["active_top5"],
            "interaction_actives_top10": i["active_top10"],
            "composite_ef1": c["ef1"], "composite_actives_top5": c["active_top5"],
            "composite_actives_top10": c["active_top10"],
            "retrieval_active_ranks": json.dumps(r["active_ranks"], sort_keys=True),
            "composite_active_ranks": json.dumps(c["active_ranks"], sort_keys=True),
            "n_composite_converged": summary["n_composite_converged"],
            "spearman_ligunity_vs_composite": summary["spearman_ligunity_vs_composite"],
        })
    frame = pd.DataFrame(rows)
    frame["delta_top5_composite_minus_retrieval"] = frame.composite_actives_top5 - frame.retrieval_actives_top5
    frame["delta_top10_composite_minus_retrieval"] = frame.composite_actives_top10 - frame.retrieval_actives_top10
    frame.to_csv(ROOT / "per_pocket_metrics.csv", index=False)

    summary = {
        "n_pockets": len(frame), "candidates_per_pocket": 285,
        "retrieved_per_pocket": 20, "poses_per_retrieved_ligand": 100,
        "total_pose_scores_expected": int(len(frame) * 20 * 100),
        "total_pose_scores_input": int(frame.n_input_poses.sum()),
        "total_pose_scores_completed": int(frame.n_scored_poses.sum()),
        "mean_full_retrieval_ef1": float(frame.full_retrieval_ef1.mean()),
        "mean_shortlist_retrieval_ef1": float(np.mean([
            json.loads((ROOT / t / "comparison_summary.json").read_text())["top20_ligunity"]["ef1"]
            for t in frame.target])),
        "mean_shortlist_interaction_ef1": float(frame.interaction_ef1.mean()),
        "mean_shortlist_composite_ef1": float(frame.composite_ef1.mean()),
        "mean_retrieval_actives_top5": float(frame.retrieval_actives_top5.mean()),
        "mean_composite_actives_top5": float(frame.composite_actives_top5.mean()),
        "mean_retrieval_actives_top10": float(frame.retrieval_actives_top10.mean()),
        "mean_composite_actives_top10": float(frame.composite_actives_top10.mean()),
        "composite_top5_improved_pockets": int((frame.delta_top5_composite_minus_retrieval > 0).sum()),
        "composite_top5_worsened_pockets": int((frame.delta_top5_composite_minus_retrieval < 0).sum()),
        "composite_top5_tied_pockets": int((frame.delta_top5_composite_minus_retrieval == 0).sum()),
        "composite_top10_improved_pockets": int((frame.delta_top10_composite_minus_retrieval > 0).sum()),
        "composite_top10_worsened_pockets": int((frame.delta_top10_composite_minus_retrieval < 0).sum()),
        "composite_top10_tied_pockets": int((frame.delta_top10_composite_minus_retrieval == 0).sum()),
        "composite_converged": int(frame.n_composite_converged.sum()),
        "composite_requested": int(len(frame) * 20),
        "selection_rule": "first 24 CASF TargetInfo targets with only public AIMNet2-supported pocket elements",
        "warning": "AIMNet2 composite is the public-checkpoint Eint + Edesolv + ELCSE reconstruction, not the unreleased trained AIMNet2(Score) checkpoint.",
    }
    (ROOT / "aggregate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
