#!/bin/bash
#SBATCH --job-name=audit-casf-t3-layers
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/xianyang/logs/audit-casf-t3-layers-%j.out

set -u
echo "Original T3 construction inputs"
for p in \
  /data/work/vs-benchmark/code/LigUnity/test_datasets/train_label_blend_seq_full.json \
  /data/work/vs-benchmark/data/raw/figshare/uniport40.clstr \
  /data/work/vs-benchmark/code/LigUnity/test_datasets/casf_label_seq.json; do
  if [ -f "$p" ]; then ls -lh "$p"; else echo "MISSING $p"; fi
done
echo "Candidate T3/CASF metadata files"
ls -lah /data/user_data/xianyang
ls -lah /data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset
ls -lah /data/user_data/xianyang/casf-2016/extracted/CASF-2016/power_scoring
ls -lah /data/user_data/xianyang/casf-2016/extracted/CASF-2016/power_screening
test -f /data/user_data/xianyang/vs-benchmark-data-2026-08-31.tar.gz && \
  tar -tzf /data/user_data/xianyang/vs-benchmark-data-2026-08-31.tar.gz | \
  grep -E '(train_label_blend_seq_full.json|uniport40.clstr|casf_label_seq.json)$' | head -n 40 || true

echo "CASF extracted root"
find /data/user_data/xianyang/casf-2016/extracted/CASF-2016 -maxdepth 3 -type f | head -n 80

echo "T3 normalized target counts"
/home/xianyang/miniconda3/envs/t3-aimnet/bin/python - <<'PY'
import pandas as pd
p = '/data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset/targets.csv.gz'
d = pd.read_csv(p)
print(d.groupby('layer').size())
print(d.columns.tolist())
print(d.head().to_string(index=False))
PY
