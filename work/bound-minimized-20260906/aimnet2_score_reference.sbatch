#!/bin/bash
#SBATCH --job-name=aimnet2-score-ref
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/xianyang/logs/aimnet2-score-reference-%j.out
#SBATCH --error=/home/xianyang/logs/aimnet2-score-reference-%j.out

set -euo pipefail

source /home/xianyang/miniconda3/etc/profile.d/conda.sh
conda activate t3-aimnet

POSES=/data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9/chunks/poses_000000_000250.sdf.gz
INTERACTION=/data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9/aimnet-smoke-000000-000250.csv
OUTPUT=/data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9/aimnet2-score-reference-000000-000250.csv
SUMMARY=/data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9/aimnet2-score-reference-000000-000250.summary.json

python /home/xianyang/aimnet2_score_pipelines/work/aimnet2_composite_smoke.py \
  --poses "$POSES" \
  --pocket /data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9/pocket.pdb \
  --interaction-csv "$INTERACTION" \
  --cpcm-model /home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt \
  --gas-model /home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt \
  --output "$OUTPUT" \
  --summary "$SUMMARY" \
  --max-poses 250 \
  --refine-steps 300 \
  --max-steps 1000 \
  --torsion-force-constant 1 \
  --device cuda

mkdir -p /home/xianyang/aimnet2_score_pipelines/work/t3-aimnet-results
cp "$OUTPUT" /home/xianyang/aimnet2_score_pipelines/work/t3-aimnet-results/
cp "$SUMMARY" /home/xianyang/aimnet2_score_pipelines/work/t3-aimnet-results/
