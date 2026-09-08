#!/bin/bash
#SBATCH --job-name=stage-casf-results
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/stage-casf-results-%j.out

set -eu

canonical=/home/xianyang/aimnet2_score_pipelines/work/casf-aimnet-multi24-results/3ebp
test -d "$canonical"
printf '%s\n' "Single-pocket staging retired; canonical 3ebp results: $canonical"
