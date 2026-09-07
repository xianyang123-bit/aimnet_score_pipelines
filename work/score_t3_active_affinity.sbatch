#!/bin/bash
#SBATCH --job-name=score-t3-active-aff
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --array=0-95%6
#SBATCH --output=/home/xianyang/logs/score-t3-active-aff-%A_%a.out

set -eu
root=/data/user_data/xianyang/t3-aimnet-affinity-balanced
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
target_dir=$($python -c "import pandas as p; print(p.read_csv('$root/tasks.tsv',sep='\\t').iloc[$SLURM_ARRAY_TASK_ID].target_dir)")

if [ ! -s "$target_dir/poses.sdf" ]; then
  echo "missing poses: $target_dir"
  exit 0
fi

$python /home/xianyang/aimnet2_score_pipelines/work/aimnet_t3_interaction.py \
  --pocket "$target_dir/pocket.pdb" --poses "$target_dir/poses.sdf" \
  --model aimnet2-2025 \
  --output "$target_dir/aimnet_interaction.csv" --batch-size 4

$python /home/xianyang/aimnet2_score_pipelines/work/aimnet2_composite_smoke.py \
  --poses "$target_dir/poses.sdf" \
  --interaction-csv "$target_dir/aimnet_interaction.csv" \
  --cpcm-model /home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt \
  --interaction-model aimnet2-2025 --gas-model /home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt \
  --output "$target_dir/aimnet2_score.csv" \
  --summary "$target_dir/aimnet2_score.summary.json" \
  --max-poses 10 --complex-max-steps 1000 --max-steps 1000 \
  --device cuda

