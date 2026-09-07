#!/bin/bash
#SBATCH --job-name=casf-ligunity
#SBATCH --partition=general
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/home/xianyang/logs/casf-ligunity-%j.out

set -eu

root=/data/user_data/xianyang/casf-2016/ligunity
envdir="$root/env"
repo=/home/xianyang/LigUnity-official
results="$root/results"
mkdir -p "$results"

export HF_HOME="$root/hf-cache"
export LIGUNITY_ESM_MODEL=facebook/esm2_t12_35M_UR50D
export PYTHONPATH="$repo"

cd "$repo"
"$envdir/bin/python" ./unimol/test.py ./test_datasets \
  --user-dir ./unimol \
  --valid-subset test \
  --results-path "$results" \
  --num-workers 4 \
  --ddp-backend=c10d \
  --batch-size 16 \
  --task test_task \
  --loss rank_softmax \
  --arch protein_ranking \
  --fp16 \
  --fp16-init-scale 4 \
  --fp16-scale-window 256 \
  --seed 1 \
  --path "$root/checkpoint/checkpoint_avg_41-50.pt" \
  --log-interval 20 \
  --log-format simple \
  --max-pocket-atoms 511 \
  --test-task CASF

"$envdir/bin/python" /home/xianyang/aimnet2_score_pipelines/work/analyze_casf_ligunity_retrieval.py \
  --results "$results" \
  --target 3ejr \
  --actives 3ejr 3dx2 3d4z 3dx1 1ps3 \
  --top-n 50
