#!/usr/bin/env bash
#SBATCH --job-name=t3-full-smina
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/xianyang/t3-full-smina-worker-%A_%a.out

set -uo pipefail

TOTAL_TASKS=487
TOTAL_WORKERS=16
failures=0

for ((task=SLURM_ARRAY_TASK_ID; task<TOTAL_TASKS; task+=TOTAL_WORKERS)); do
  echo "[$(date -Is)] worker=${SLURM_ARRAY_TASK_ID} task=${task}"
  if ! /home/xianyang/miniconda3/envs/t3-aimnet/bin/python \
      /home/xianyang/aimnet2_score_pipelines/work/t3_full_dock.py \
      --data /data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset \
      --work /data/user_data/xianyang/t3-aimnet-full/full-dock-L4 \
      --env-bin /home/xianyang/miniconda3/envs/t3-aimnet/bin \
      dock-task --task-id "${task}" --cpus "${SLURM_CPUS_PER_TASK}"; then
    failures=$((failures + 1))
    echo "[$(date -Is)] task=${task} FAILED; continuing"
  fi
done

echo "[$(date -Is)] worker=${SLURM_ARRAY_TASK_ID} finished failures=${failures}"
exit "${failures}"
