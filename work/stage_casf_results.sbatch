#!/bin/bash
#SBATCH --job-name=stage-casf-results
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/stage-casf-results-%j.out

set -eu
src=/data/user_data/xianyang/casf-2016/ligunity
rerank="$src/rerank-3ebp-top20"
out=/home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results
mkdir -p "$out"
cp "$src/results/3ebp_ligunity_all285.csv" "$out/"
cp "$src/results/3ebp_ligunity_top50.csv" "$out/"
cp "$src/results/casf_target_support_and_retrieval.csv" "$out/"
cp "$src/results/selected_supported_target.json" "$out/"
cp "$rerank/preparation.json" "$out/"
cp "$rerank/aimnet_interaction_allposes.csv" "$out/"
cp "$rerank/aimnet_interaction_bestpose.csv" "$out/"
cp "$rerank/aimnet_interaction_bestpose.summary.json" "$out/"
cp "$rerank/aimnet_interaction_bestpose.sdf" "$out/"
cp "$rerank/aimnet2_composite_rerank.csv" "$out/"
cp "$rerank/aimnet2_composite_rerank.summary.json" "$out/"
cp "$rerank/ranking_comparison_all20.csv" "$out/"
cp "$rerank/comparison_summary.json" "$out/"
cp /home/xianyang/logs/casf-ligunity-10318400.out "$out/ligunity-run.log"
cp /home/xianyang/logs/casf-aimnet-rerank-10318481.out "$out/aimnet-rerank-run.log"
du -sh "$out"
find "$out" -maxdepth 1 -type f -printf '%s %f\n' | sort -n
