#!/bin/bash
#SBATCH --job-name=prep-t3-active-aff
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/xianyang/logs/prep-t3-active-aff-%j.out

set -eu
: "${PREPARED_POCKET_ROOT:?Set PREPARED_POCKET_ROOT to the validated LAYER/UNIPROT receptor tree}"
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/prepare_t3_active_affinity.py \
  --prepared-pocket-root "$PREPARED_POCKET_ROOT" \
  --data /data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset \
  --work "${PREPARED_WORK_ROOT:?Set PREPARED_WORK_ROOT to a new output directory}" \
  --env-bin /home/xianyang/miniconda3/envs/t3-aimnet/bin \
  --targets-per-layer 24 --ligands-per-target 10 --seed 20260904

