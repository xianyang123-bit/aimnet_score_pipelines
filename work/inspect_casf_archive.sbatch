#!/bin/bash
#SBATCH --job-name=inspect-casf
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --output=/home/xianyang/logs/inspect-casf-%j.out
#SBATCH --error=/home/xianyang/logs/inspect-casf-%j.out

set -euo pipefail
python3 /home/xianyang/aimnet2_score_pipelines/work/inspect_casf_archive.py \
  /data/user_data/xianyang/casf-2016/casf_2016.tar.gz
