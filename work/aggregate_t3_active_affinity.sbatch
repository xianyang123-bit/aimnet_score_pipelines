#!/bin/bash
#SBATCH --job-name=agg-t3-active-aff
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/home/xianyang/logs/agg-t3-active-aff-%j.out

set -eu
root=/data/user_data/xianyang/t3-aimnet-affinity-balanced
out=/home/xianyang/aimnet2_score_pipelines/work/t3-aimnet-affinity-balanced-results
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
$python /home/xianyang/aimnet2_score_pipelines/work/aggregate_t3_active_affinity.py --work "$root"
mkdir -p "$out"
cp "$root/manifest.json" "$root/selected_targets.csv" "$root/tasks.tsv" \
   "$root/per_target_metrics.csv" "$root/per_ligand_scores.csv.gz" \
   "$root/aggregate_summary.json" "$out/"
find "$root" -name dock_status.json -print0 | xargs -0 -r cp --parents -t "$out"
du -sh "$out"
