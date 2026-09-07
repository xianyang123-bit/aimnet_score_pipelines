# AIMNet2 score pipelines

This is the canonical Babel package for the scripts and outputs created during
the T3/CASF docking, retrieval, and AIMNet2 scoring work. Its deployed path is:

```text
/home/xianyang/aimnet2_score_pipelines
```

Package-internal paths in the Python and SLURM files point to
`/home/xianyang/aimnet2_score_pipelines/work`. Large datasets, environments,
model checkpoints, and repositories remain at their external Babel paths.

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
- `work/casf-aimnet-results/` and `work/casf-aimnet-multi24-results/`: staged
  CASF results.
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

## Numerical correction, 2026-09-06

The T3 affinity, CASF reranking, and 250-ligand smoke-test results were recalculated using explicit reference neighbor lists and full-range Coulomb neighbors for the legacy checkpoint. Scorers now import work/aimnet_safe_calculator.py. See work/correction-20260906/report.md for metrics, validation, limitations, and original-result backups. Eight historical exploratory CSVs were not regenerated; see legacy_experiments_not_recomputed.json. This correction retains the original checkpoints and does not turn the reconstructed score into an official AIMNet2(Score) release.

Outdated result backups and historical result CSVs were deleted on 2026-09-06 at the user's request. Current results and comparison audits remain; see work/correction-20260906/cleanup_log.json.

Obsolete code backups and one-time correction/debug scripts were deleted at the user's request. Current pipeline and rerun utilities remain; see work/correction-20260906/code_cleanup_log.json.
