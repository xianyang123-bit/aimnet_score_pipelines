#!/bin/bash
#SBATCH --job-name=inspect-casf-screen
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/inspect-casf-screen-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016/extracted/CASF-2016
archive=/data/user_data/xianyang/casf-2016/casf_2016.tar.gz

echo TARGET_INFO
sed -n '1,120p' "$root/power_screening/TargetInfo.dat"

echo LIGAND_INFO
sed -n '1,45p' "$root/power_screening/LigandInfo.dat"

echo SCREENING_TOP_LEVEL
/home/xianyang/miniconda3/bin/python /home/xianyang/aimnet2_score_pipelines/work/inspect_casf_screening_archive.py "$archive"
