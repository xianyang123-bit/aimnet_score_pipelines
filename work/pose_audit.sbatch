#!/bin/bash
#SBATCH --job-name=pose-audit
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home/xianyang/logs/pose-audit-%j.out
#SBATCH --error=/home/xianyang/logs/pose-audit-%j.out

set -euo pipefail

ORIGINAL=/data/yicheng/xqc/vs-benchmark/dock

echo ORIGINAL_T6_DOCK_DIRECTORY
ls -ld "$ORIGINAL" 2>&1 || true

if [[ -d "$ORIGINAL" ]]; then
  echo ORIGINAL_T6_POSE_FILES
  find "$ORIGINAL" -type f \( \
    -iname 'poses.sdf' -o -iname 'poses.sdf.gz' -o \
    -iname '*docked*.sdf' -o -iname '*docked*.sdf.gz' \
  \) -printf '%p|%s bytes\n'

  echo ORIGINAL_T6_TARGET_DIRECTORIES
  find "$ORIGINAL" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
fi
