#!/bin/bash
#SBATCH --job-name=analyze-casf-rerank
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/analyze-casf-rerank-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016/ligunity/rerank-3ebp-top20
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/analyze_casf_aimnet_rerank.py
# Analysis remains under "$root"; canonical staged results are in
# /home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-multi24-results/3ebp
