#!/bin/bash
#SBATCH --job-name=casf-t3-overlap
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home/xianyang/logs/casf-t3-overlap-%j.out

set -eu
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python /home/xianyang/aimnet2_score_pipelines/work/audit_casf_t3_overlap.py
