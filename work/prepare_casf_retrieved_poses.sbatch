#!/bin/bash
#SBATCH --job-name=prep-casf-rerank
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/xianyang/logs/prep-casf-rerank-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016
out="$root/ligunity/rerank-3ebp-top20"
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python \
  /home/xianyang/aimnet2_score_pipelines/work/prepare_casf_retrieved_poses.py \
  --top-csv "$root/ligunity/results/3ebp_ligunity_top50.csv" \
  --archive "$root/casf_2016.tar.gz" \
  --casf-root "$root/extracted/CASF-2016" \
  --target 3ebp \
  --top-n 20 \
  --output "$out" \
  --obabel /home/xianyang/miniconda3/envs/obabel/bin/obabel
