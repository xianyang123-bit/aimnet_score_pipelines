#!/bin/bash
#SBATCH --job-name=audit-ligunity-casf
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home/xianyang/logs/audit-ligunity-casf-%j.out

set -u

echo ROOTS
ls -ld \
  /data/user_data/xianyang/t3-aimnet-full/vs-benchmark \
  /data/user_data/xianyang/casf-2016/extracted/CASF-2016

echo ASSETS
find \
  /data/user_data/xianyang/t3-aimnet-full/vs-benchmark \
  -maxdepth 7 -type f \
  \( -iname '*.pt' -o -iname '*.pth' -o -iname '*.lmdb' -o -iname '*casf*' \) \
  -printf '%s %p\n' 2>/dev/null | sort -nr | head -250

echo CASF_LAYOUT
find /data/user_data/xianyang/casf-2016/extracted/CASF-2016 \
  -maxdepth 3 -type f -printf '%p\n' | head -120
