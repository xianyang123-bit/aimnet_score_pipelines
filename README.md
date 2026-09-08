# AIMNet2 score pipelines

This is the canonical Babel package for the scripts and outputs created during
the T3/CASF docking, retrieval, and AIMNet2 scoring work. Its deployed path is:

```text
/home/xianyang/aimnet2_score_pipelines
```

Package-internal paths in the Python and SLURM files point to
`/home/xianyang/aimnet2_score_pipelines/work`. Large datasets, environments,
model checkpoints, and repositories remain at their external Babel paths.

<!-- AIMNET2025_RESULTS_START -->
## Latest benchmark results: AIMNet2(2025) member 0 (2026-09-06)

These are the current results, superseding the earlier checkpoint tables. Interaction energies and ligand minimization within a fixed protein pocket use AIMNet2(2025) member 0 (aimnet2-b973c-2025-d3_0). Desolvation retains the matched legacy wB97M gas / wB97M-CPCM pair. Composite = minimized interaction + desolvation + local ligand strain; lower scores are better.

### T3 affinity ranking

917 ligand cases across 93 target/layer systems were scored. Values below are the mean per-system Spearman correlation between negative energy and pAffinity; larger is better. SMINA uses negative docking affinity. This is an active-only affinity-ranking benchmark.

| Layer | Systems | Ligands | Interaction before minimization | Interaction after minimization | Composite | SMINA | Previous composite |
|---|---:|---:|---:|---:|---:|---:|---:|
| L1 | 24 | 231 | -0.0534 | -0.0407 | -0.0877 | 0.0302 | -0.0165 |
| L2 | 24 | 236 | 0.1070 | 0.0478 | 0.1018 | 0.0945 | 0.0441 |
| L3 | 23 | 230 | 0.0854 | 0.2350 | 0.2174 | 0.1069 | 0.0698 |
| L4 | 22 | 220 | -0.0657 | 0.1330 | 0.1169 | -0.0104 | 0.0177 |

The composite improves over the previous checkpoint in L2–L4 and declines in L1. L4 has 21 evaluable correlations from 22 scored systems because one system has constant affinity. Three empty pose inputs are excluded: L3/P18564, L4/Q15726, and L4/Q8TDU9.

**T3 retrieval comparison:** matched retrieval scores are absent from these scored T3 outputs, so a T3 retrieval Spearman baseline cannot be reported from this bundle. SMINA is a docking baseline, not the retrieval model. The direct retrieval comparison below uses the CASF shortlist and a different metric.

### CASF-2016 subset: comparison with LigUnity retrieval

All 47,896 available poses for 480 candidate ligand cases across 24 pockets were rescored. The lowest interaction-energy pose per ligand was selected and minimized; selections changed for 320/480 cases. The separate 3ebp example (2,000 poses, 20 ligand cases) is excluded from these 24-pocket means.

| Ranking method | Mean shortlist EF1 | Mean actives in top 5 | Mean actives in top 10 |
|---|---:|---:|---:|
| LigUnity retrieval | 3.2024 | 4.3333 | 4.7917 |
| AIMNet2(2025) minimized interaction | 0.9722 | 1.2917 | 2.7083 |
| AIMNet2(2025) composite | 0.9722 | 1.4167 | 2.4583 |

Interaction ranking before minimization gives mean shortlist EF1 = 1.6389. Minimization reduces enrichment in this subset. Retrieval remains stronger than either minimized AIMNet2 ranking.

Each pocket contains 20 candidates selected by retrieval. EF1 uses the top ceil(0.01 × 20) = 1 candidate and is normalized by that shortlist's active fraction. These are conditional reranking results on retrieval-selected shortlists, not full-library CASF screening metrics. Previous minimized interaction/composite EF1 values were 0.7579/0.9722; the new values are 0.9722/0.9722. The comparison with the previous checkpoint also includes changed pose selection.

### T3 smoke set

250 nonrandomly selected ligands (58 actives, 192 decoys). These AUROC values measure screening discrimination, not the active-only affinity correlation above.

| Method | AUROC | EF1 |
|---|---:|---:|
| Interaction before minimization | 0.5796 | 0.0000 |
| Minimized interaction | 0.5422 | 1.4368 |
| Composite | 0.5216 | 0.0000 |

Previous minimized interaction/composite AUROC values were 0.4543/0.4472. Both improve, but remain close to chance. Matched retrieval scores are absent from the smoke inputs, so no smoke retrieval AUROC is claimed.

### Validation and interpretation

All 1,667 ligand cases scored without errors, input IDs matched, and component-energy arithmetic passed at 1e-6 kcal/mol. Fixed-pocket minimization converged for 1,440 cases; isolated CPCM minimization for 1,490; both for 1,313. Finite nonconverged results are retained and flagged. No final complex energy exceeded its initial value.

This is a reconstruction of the composite energy expression using public checkpoints, not an official affinity-trained AIMNet2(Score) release. Improved reference interaction-energy accuracy does not guarantee improved affinity or screening ranking.

Full tables, minimized poses, scores, and provenance: [current results](work/aimnet2025-benchmarks-20260906/README.md).
<!-- AIMNET2025_RESULTS_END -->

## Contents

- `work/prepare_t3_active_affinity.py` through
  `work/aggregate_t3_active_affinity.py`: balanced T3 L1-L4 affinity-ranking
  pipeline.
- `work/*t3_active_affinity*.sbatch`: preparation, docking, GPU scoring,
  aggregation, and artifact-staging jobs.
- `work/t3-aimnet-affinity-balanced-results/`: compact T3 tables and metadata.
- `work/t3-aimnet-affinity-balanced-poses-and-scores.tar.gz`: per-target
  pockets, ligand inputs, docked `poses.sdf` files, and AIMNet2 component scores.
- `work/*casf*` and `work/*ligunity*`: CASF-2016 retrieval, pose preparation,
  reranking, multi-pocket evaluation, and audit scripts.
- `work/casf-aimnet-multi24-results/`: canonical staged CASF results,
  including the original `3ebp` single-pocket case.
- `work/t3_full_dock*`, `work/aimnet_*chunk*`, and
  `work/aimnet2_score_reference*`: earlier T3 docking and scoring experiments.
- `MANIFEST.sha256`: SHA-256 checksums for every file in this consolidated copy.

## External data and models not duplicated

The large datasets, model checkpoints, environments, and repositories remain on
Babel. Existing scripts refer to these paths directly:

```text
/data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset
/data/user_data/xianyang/t3-aimnet-affinity-balanced
/data/user_data/xianyang/casf-2016
/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt
/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt
/home/xianyang/miniconda3/envs/t3-aimnet
/home/xianyang/miniconda3/envs/obabel
/home/xianyang/LigUnity-official
```

The original uploaded benchmark archive is retained at:

```text
/data/user_data/xianyang/vs-benchmark-data-2026-08-31.tar.gz
```

## Important scoring caveat

The included composite scorer reconstructs the published AIMNet2(Score) energy
expression as interaction energy + desolvation + ligand conformational strain
using public AIMNet2 checkpoints. It is not an official released
AIMNet2(Score) implementation or affinity-trained checkpoint.

The SLURM files contain absolute Babel paths. Update the package root if this
directory is deployed under a different account or location.

## Historical numerical correction, 2026-09-06

The T3 affinity, CASF reranking, and 250-ligand smoke-test results were recalculated using explicit reference neighbor lists and full-range Coulomb neighbors for the legacy checkpoint. Scorers now import work/aimnet_safe_calculator.py. See work/correction-20260906/report.md for metrics, validation, limitations, and the retained correction audit record. Eight historical exploratory CSVs were not regenerated; see legacy_experiments_not_recomputed.json. This correction retains the original checkpoints and does not turn the reconstructed score into an official AIMNet2(Score) release.

Outdated result backups and historical result CSVs were deleted on 2026-09-06 at the user's request. Current results and comparison audits remain; see work/correction-20260906/cleanup_log.json.

Obsolete code backups and one-time correction/debug scripts were deleted at the user's request. Current pipeline and rerun utilities remain; see work/correction-20260906/code_cleanup_log.json.
