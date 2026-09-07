#!/usr/bin/env bash
#SBATCH --job-name=t3-full-smina
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/home/xianyang/t3-full-smina-%A_%a.out

set -euo pipefail

/home/xianyang/miniconda3/envs/t3-aimnet/bin/python \
  /home/xianyang/aimnet2_score_pipelines/work/t3_full_dock.py \
  --data /data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset \
  --work /data/user_data/xianyang/t3-aimnet-full/full-dock-L4 \
  --env-bin /home/xianyang/miniconda3/envs/t3-aimnet/bin \
  dock-task --cpus "${SLURM_CPUS_PER_TASK}"
