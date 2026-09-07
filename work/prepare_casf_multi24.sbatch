#!/bin/bash
#SBATCH --job-name=prep-casf-multi24
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/xianyang/logs/prep-casf-multi24-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/prepare_casf_multi24.py \
  --casf-root "$root/extracted/CASF-2016" \
  --archive "$root/casf_2016.tar.gz" \
  --retrieval-results "$root/ligunity/results" \
  --support-csv "$root/ligunity/results/casf_target_support_and_retrieval.csv" \
  --output "$root/ligunity/multi24" \
  --obabel /home/xianyang/miniconda3/envs/obabel/bin/obabel \
  --n-targets 24 \
  --top-n 20
