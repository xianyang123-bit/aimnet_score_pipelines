#!/bin/bash
#SBATCH --job-name=select-casf-target
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/select-casf-target-%j.out

set -eu
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python \
  /home/xianyang/aimnet2_score_pipelines/work/select_supported_casf_target.py \
  --casf-root /data/user_data/xianyang/casf-2016/extracted/CASF-2016 \
  --results /data/user_data/xianyang/casf-2016/ligunity/results
