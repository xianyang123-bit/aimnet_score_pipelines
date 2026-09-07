#!/bin/bash
#SBATCH --job-name=casf-aimnet-rerank
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home/xianyang/logs/casf-aimnet-rerank-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016/ligunity/rerank-3ebp-top20
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
poses="$root/3ebp_retrieved_top20_allposes.sdf"

"$python" /home/xianyang/aimnet2_score_pipelines/work/aimnet_casf_interaction.py \
  --pocket "$root/pocket.pdb" \
  --poses "$poses" \
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

mkdir -p /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results
cp "$root/aimnet_interaction_bestpose.csv" /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
cp "$root/aimnet_interaction_bestpose.summary.json" /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
cp "$root/aimnet2_composite_rerank.csv" /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
cp "$root/aimnet2_composite_rerank.summary.json" /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
cp /data/user_data/xianyang/casf-2016/ligunity/results/3ebp_ligunity_all285.csv /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
cp /data/user_data/xianyang/casf-2016/ligunity/results/selected_supported_target.json /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-results/
