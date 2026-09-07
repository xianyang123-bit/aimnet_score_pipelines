#!/bin/bash
#SBATCH --job-name=aggregate-casf-multi24
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home/xianyang/logs/aggregate-casf-multi24-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016/ligunity/multi24
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
"$python" /home/xianyang/aimnet2_score_pipelines/work/aggregate_casf_multi24.py
out=/home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-multi24-results
mkdir -p "$out"
cp "$root/selected_targets.csv" "$out/"
cp "$root/per_pocket_metrics.csv" "$out/"
cp "$root/aggregate_summary.json" "$out/"
for target_dir in "$root"/*; do
  [ -d "$target_dir" ] || continue
  target=$(basename "$target_dir")
  mkdir -p "$out/$target"
  cp "$target_dir/comparison_summary.json" "$out/$target/"
  cp "$target_dir/ranking_comparison_all20.csv" "$out/$target/"
  cp "$target_dir/aimnet2_composite_rerank.csv" "$out/$target/"
  cp "$target_dir/aimnet2_composite_rerank.summary.json" "$out/$target/"
  cp "$target_dir/aimnet_interaction_bestpose.csv" "$out/$target/"
  cp "$target_dir/aimnet_interaction_bestpose.summary.json" "$out/$target/"
done
du -sh "$out"
