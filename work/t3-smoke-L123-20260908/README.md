# T3 L1-L3 stratified smoke screening

Scope: up to 250 ligands per target, without replacement, using both original active and decoy labels. All targets in T3 L1-L3 are inventoried. Unavailable pockets are explicitly reported. Samples preserve source prevalence as closely as integer counts permit, typically 5 actives and 245 decoys; this is not the previous active-enriched first chunk.

Reproducibility: base seed 20260908; each target receives the first 8 SHA256 bytes of seed:layer:uniprot as a little-endian integer. Class pools are sorted by mol_id, sampled using NumPy default_rng, and the combined sample is shuffled. library.csv.gz retains original source indices, labels and IDs. The source full-library preparation is read-only.

Protocol: frozen code from aimnet2025-benchmarks-20260906 via the cancelled full-run snapshot. SMINA exhaustiveness 8, one pose, seed 1; 6A pocket box plus 4A padding. AIMNet2(2025) member 0 fixed-pocket minimization, 1000 steps, 0.002 eV/A; legacy wB97M gas/CPCM desolvation and local strain, 1000 steps. Same reconstructed scoring protocol as the earlier smoke result. Original heavy-atom pocket and charge limitations remain.

Metrics: SMINA, pre-minimization interaction, minimized interaction and composite. Lower energy is better. AUROC uses average ranks. EF@1% uses k=ceil(0.01*N), so k=3 for 250 ligands; expected active hits across boundary ties prevent input-order bias. EF is fold enrichment. N and A refer to the sampled library, not the full source library.

All sampled IDs remain in primary sample_failure_last metrics; failed scores are tied last. scored_only metrics and per-method coverage are also reported. All-failed targets receive NA metrics. Finite nonconverged scores are retained and flagged. Layer means are unweighted macro means over completed evaluable targets, with pending, failed and unavailable counts. No cross-target pooling of energies.

Execution: preparation validates class counts, label alignment and sampling reproducibility, then submits a two-ligand GPU pilot. Eight GPU workers start only after it passes. Checkpoints are every 25 ligands, with dependent continuation jobs before the wall-time limit. Target failures are recorded and other targets continue; repeated failures stop the lane for inspection.

Summaries in this directory: manifest.json, target_inventory.csv, per_target_metrics.csv, per_layer_metrics.csv, summary.json, pilot.json, launch.json, worker_*.json. Large poses and scores: /data/user_data/xianyang/t3-smoke-L123-20260908 (compute nodes). The full-library run was cancelled at user request before full workers started.

## Faster execution

fast.py replaces the sequential workers with 16 four-CPU docking workers and 8 GPU allocations. A four-ligand same-GPU serial/concurrent benchmark selects one or two scoring processes per GPU and checks energy/flag agreement. Scientific model code, sampling, SMINA settings and minimization limits are unchanged. babel-x5-20 is excluded after CUDA initialization and memory failures. Completed chunks and successfully docked poses are reused; interrupted unfinished chunks can restart. Tasks interleave L1/L2/L3 targets. CPU and GPU continuations checkpoint separately. Monitor fast/jobs.json, fast/cpu_*.json, fast/score_*.json, fast/gpu_*.json, fast/errors and fast/logs. Original status files are archived in fast/legacy_status; summary CSV paths stay unchanged.

Concurrent scoring was rejected by the numerical-equivalence check (maximum component difference about 0.295 kcal/mol). Production uses one scorer per GPU, with the original four CPU threads and 25-ligand scoring order. Only docking/scoring stage separation is enabled. The serial four-ligand reference completed successfully; the concurrent timings are not an approved speedup.
