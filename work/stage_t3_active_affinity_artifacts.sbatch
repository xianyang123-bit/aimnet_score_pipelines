#!/bin/bash
#SBATCH --job-name=stage-t3-aff-artifacts
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/home/xianyang/logs/stage-t3-aff-artifacts-%j.out

set -eu
root=/data/user_data/xianyang/t3-aimnet-affinity-balanced
archive=/home/xianyang/aimnet2_score_pipelines/work/t3-aimnet-affinity-balanced-poses-and-scores.tar.gz
cd "$root"
{
  printf '%s\n' manifest.json selected_targets.csv tasks.tsv \
    per_target_metrics.csv per_ligand_scores.csv.gz aggregate_summary.json
  find L1 L2 L3 L4 -type f \( \
    -name pocket.pdb -o -name pocket.pdbqt -o -name box.json -o \
    -name ligands.sdf -o -name poses.sdf -o -name dock_status.json -o \
    -name smina.log -o -name aimnet_interaction.csv -o \
    -name aimnet2_score.csv -o -name aimnet2_score.summary.json \
  \) -print
} | tar -czf "$archive" -T -
ls -lh "$archive"
