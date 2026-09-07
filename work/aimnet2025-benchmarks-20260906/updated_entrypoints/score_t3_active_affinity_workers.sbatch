#!/bin/bash
#SBATCH --job-name=score-t3-aff-worker
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=/home/xianyang/logs/score-t3-aff-worker-%A_%a.out

set -u
root=/data/user_data/xianyang/t3-aimnet-affinity-balanced
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
total=96
workers=6
failures=0
for ((task=SLURM_ARRAY_TASK_ID; task<total; task+=workers)); do
  target_dir=$($python -c "import pandas as p; print(p.read_csv('$root/tasks.tsv',sep='\\t').iloc[$task].target_dir)")
  if [ ! -s "$target_dir/poses.sdf" ]; then
    echo "missing poses task=$task path=$target_dir"
    failures=$((failures+1))
    continue
  fi
  echo "[$(date -Is)] score task=$task path=$target_dir"
  if ! $python /home/xianyang/aimnet2_score_pipelines/work/aimnet_t3_interaction.py \
      --pocket "$target_dir/pocket.pdb" --poses "$target_dir/poses.sdf" \
      --model aimnet2-2025 \
      --output "$target_dir/aimnet_interaction.csv" --batch-size 4; then
    failures=$((failures+1)); continue
  fi
  if ! $python /home/xianyang/aimnet2_score_pipelines/work/aimnet2_composite_smoke.py \
      --poses "$target_dir/poses.sdf" --interaction-csv "$target_dir/aimnet_interaction.csv" \
      --cpcm-model /home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt \
      --interaction-model aimnet2-2025 --gas-model /home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt \
      --output "$target_dir/aimnet2_score.csv" --summary "$target_dir/aimnet2_score.summary.json" \
      --max-poses 10 --complex-max-steps 1000 --max-steps 1000 \
      --device cuda; then
    failures=$((failures+1))
  fi
done
echo "worker=${SLURM_ARRAY_TASK_ID} failures=$failures"
exit 0
