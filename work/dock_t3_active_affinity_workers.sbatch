#!/bin/bash
#SBATCH --job-name=dock-t3-aff-worker
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --array=0-15
#SBATCH --output=/home/xianyang/logs/dock-t3-aff-worker-%A_%a.out

set -u
total=96
workers=16
failures=0
for ((task=SLURM_ARRAY_TASK_ID; task<total; task+=workers)); do
  echo "[$(date -Is)] dock task=$task"
  /home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/dock_t3_active_task.py \
    --work /data/user_data/xianyang/t3-aimnet-affinity-balanced \
    --env-bin /home/xianyang/miniconda3/envs/t3-aimnet/bin \
    --task-id "$task" --cpus "${SLURM_CPUS_PER_TASK}" --exhaustiveness 8 || failures=$((failures+1))
done
echo "worker=${SLURM_ARRAY_TASK_ID} failures=$failures"
exit 0

