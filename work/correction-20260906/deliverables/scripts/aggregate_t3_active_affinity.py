#!/usr/bin/env python3
"""Aggregate active-only T3 AIMNet2(Score) affinity-ranking results by layer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr


def corr(x, y, fn):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(fn(x[ok], y[ok]).statistic) if ok.sum() >= 3 and np.std(x[ok]) > 0 and np.std(y[ok]) > 0 else float("nan")


def bootstrap_ci(values, seed=20260904, n_boot=2000):
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    if not len(v): return [None, None]
    rng = np.random.default_rng(seed)
    means = np.mean(rng.choice(v, (n_boot, len(v)), replace=True), axis=1)
    return [float(np.quantile(means, .025)), float(np.quantile(means, .975))]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    args = ap.parse_args()
    root = Path(args.work)
    tasks = pd.read_csv(root / "tasks.tsv", sep="\t")
    ligand_frames, target_rows = [], []
    for task in tasks.itertuples(index=False):
        d = Path(task.target_dir)
        score_path = d / "aimnet2_score.csv"
        if not score_path.exists():
            target_rows.append({"layer": task.layer, "uniprot": task.uniprot, "n_scored": 0})
            continue
        frame = pd.read_csv(score_path)
        frame["layer"] = task.layer; frame["uniprot"] = task.uniprot
        valid = frame[frame.get("error", "").fillna("") == ""].copy() if "error" in frame else frame.copy()
        valid = valid[np.isfinite(valid.aimnet2_composite_kcal_mol) & np.isfinite(valid.paff)]
        ligand_frames.append(valid)
        aim_score = -valid.aimnet2_composite_kcal_mol
        int_score = -valid.aimnet_interaction_kcal_mol
        smina_score = valid.smina_score
        target_rows.append({
            "layer": task.layer, "uniprot": task.uniprot, "n_scored": len(valid),
            "aimnet_spearman": corr(aim_score, valid.paff, spearmanr),
            "aimnet_kendall": corr(aim_score, valid.paff, kendalltau),
            "aimnet_pearson": corr(aim_score, valid.paff, pearsonr),
            "interaction_spearman": corr(int_score, valid.paff, spearmanr),
            "smina_spearman": corr(smina_score, valid.paff, spearmanr),
            "paff_range": float(valid.paff.max() - valid.paff.min()) if len(valid) else float("nan"),
            "optimization_converged": int(valid.optimization_converged.fillna(False).sum()) if "optimization_converged" in valid else 0,
        })
    per_target = pd.DataFrame(target_rows)
    per_target.to_csv(root / "per_target_metrics.csv", index=False)
    all_ligands = pd.concat(ligand_frames, ignore_index=True) if ligand_frames else pd.DataFrame()
    all_ligands.to_csv(root / "per_ligand_scores.csv.gz", index=False)
    layers = {}
    for layer in ("L1", "L2", "L3", "L4"):
        f = per_target[(per_target.layer == layer) & per_target.aimnet_spearman.notna()]
        layers[layer] = {
            "n_targets": len(f), "n_ligands": int(f.n_scored.sum()),
            "mean_aimnet_spearman": float(f.aimnet_spearman.mean()),
            "median_aimnet_spearman": float(f.aimnet_spearman.median()),
            "aimnet_spearman_ci95": bootstrap_ci(f.aimnet_spearman),
            "mean_aimnet_kendall": float(f.aimnet_kendall.mean()),
            "mean_aimnet_pearson": float(f.aimnet_pearson.mean()),
            "fraction_aimnet_spearman_positive": float((f.aimnet_spearman > 0).mean()),
            "mean_interaction_spearman": float(f.interaction_spearman.mean()),
            "mean_smina_spearman": float(f.smina_spearman.mean()),
        }
    l1, l4 = layers["L1"]["mean_aimnet_spearman"], layers["L4"]["mean_aimnet_spearman"]
    summary = {
        "scope": "balanced active-only T3 pilot; no decoys, retrieval, reranking, or EF metrics",
        "score_definition": "AIMNet2(Score) public-checkpoint reconstruction: Eint + Edesolv + ELCSE; lower is stronger",
        "layers": layers,
        "l1_to_l4_absolute_change": l4 - l1,
        "l1_to_l4_relative_change": (l4 - l1) / abs(l1) if l1 else None,
        "warnings": [
            "This is not an official released AIMNet2(Score) implementation/checkpoint.",
            "T3 LMDB pockets omit residue identities and hydrogens; pocket charge is assumed zero.",
            "Ligand poses are single smina poses, so pose error can dominate affinity ranking.",
        ],
    }
    (root / "aggregate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(per_target.to_string(index=False))


if __name__ == "__main__":
    main()
