#!/bin/bash
#SBATCH --job-name=casf-multi24-aimnet
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --array=0-23%6
#SBATCH --output=/home/xianyang/logs/casf-multi24-aimnet-%A_%a.out

set -eu
multi=/data/user_data/xianyang/casf-2016/ligunity/multi24
target=$(/home/xianyang/miniconda3/envs/t3-aimnet/bin/python -c "import pandas as p; print(p.read_csv('$multi/selected_targets.csv').iloc[$SLURM_ARRAY_TASK_ID].target)")
root="$multi/$target"
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python

"$python" /home/xianyang/aimnet2_score_pipelines/work/aimnet_casf_interaction.py \
  --pocket "$root/pocket.pdb" \
  --poses "$root/${target}_top20_allposes.sdf" \
  --model aimnet2-2025 \
  --output-all "$root/aimnet_interaction_allposes.csv" \
  --output-best "$root/aimnet_interaction_bestpose.csv" \
  --output-best-sdf "$root/aimnet_interaction_bestpose.sdf" \
  --batch-size 4

"$python" /home/xianyang/aimnet2_score_pipelines/work/aimnet2_composite_smoke.py \
  --poses "$root/aimnet_interaction_bestpose.sdf" \
  --interaction-csv "$root/aimnet_interaction_bestpose.csv" \
  --cpcm-model /home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt \
  --interaction-model aimnet2-2025 --gas-model /home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt \
  --output "$root/aimnet2_composite_rerank.csv" \
  --summary "$root/aimnet2_composite_rerank.summary.json" \
  --max-poses 20 \
  --complex-max-steps 1000 \
  --max-steps 1000 \
  --device cuda

"$python" /home/xianyang/aimnet2_score_pipelines/work/analyze_casf_target.py --target-dir "$root" --target "$target"
