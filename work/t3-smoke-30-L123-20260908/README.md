# T3 30-ligand smoke run
Exactly 5 actives and 25 decoys are sampled without replacement per L1-L3 target from the preceding reproducible 250-ligand sample. Seed is derived from 20260908:30:layer:uniprot. Sampling never uses prediction scores or job completion. Existing validated scores and completed docking poses are reused by ligand ID.

The source snapshot and original scoring protocol are retained: SMINA exhaustiveness 8, seed 1, one pose; AIMNet2(2025) member 0 fixed-pocket minimization; original gas/CPCM terms and 1000-step limits. One scoring process per GPU with four CPU threads. No new concurrency benchmark is launched. Checkpoints are ten ligands.

Primary EF@1% and AUROC use all 30 sampled IDs with failed scores tied last; scored-only metrics and coverage are also reported. EF@1% uses the top one ligand and baseline active prevalence 5/30. Without boundary ties, EF@1% is 0 or 6. Layer means exclude unavailable/pending targets and retain their counts.

Execution uses 16 CPU docking workers and eight GPU workers, excluding babel-x5-20 after earlier CUDA failures. Jobs and progress are under fast/, with per_target_metrics.csv, per_layer_metrics.csv and summary.json in this directory. Large data are under /data/user_data/xianyang/t3-smoke-30-L123-20260908 on compute nodes. The preceding 250-ligand jobs have been cancelled.
