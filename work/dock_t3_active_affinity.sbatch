#!/bin/bash
#SBATCH --job-name=dock-t3-active-aff
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --array=0-95%24
#SBATCH --output=/home/xianyang/logs/dock-t3-active-aff-%A_%a.out

set -eu
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/dock_t3_active_task.py \
  --work /data/user_data/xianyang/t3-aimnet-affinity-balanced \
  --env-bin /home/xianyang/miniconda3/envs/t3-aimnet/bin \
  --cpus "${SLURM_CPUS_PER_TASK}" --exhaustiveness 8

