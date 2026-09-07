#!/bin/bash
#SBATCH --job-name=stage-casf-example
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/stage-casf-example-%j.out

set -eu
root=/data/user_data/xianyang/casf-2016
out="$root/ligunity/example"
mkdir -p "$out"
cp "$root/extracted/CASF-2016/README.md" /home/xianyang/aimnet2_score_pipelines/work/CASF-2016-README.md
tar -xzf "$root/casf_2016.tar.gz" -C "$out" \
  CASF-2016/decoys_screening/3ejr/3ejr_3ejr.mol2
ls -lh "$out/CASF-2016/decoys_screening/3ejr/3ejr_3ejr.mol2"
sed -n '1,100p' "$out/CASF-2016/decoys_screening/3ejr/3ejr_3ejr.mol2"
