#!/bin/bash
#SBATCH --job-name=extract-casf
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/home/xianyang/logs/extract-casf-%j.out
#SBATCH --error=/home/xianyang/logs/extract-casf-%j.out

set -euo pipefail

ROOT=/data/user_data/xianyang/casf-2016
ARCHIVE="$ROOT/casf_2016.tar.gz"
DEST="$ROOT/extracted"
EXPECTED_MD5=b433ef3da03ae94731d9bd0b67df8138

echo "$EXPECTED_MD5  $ARCHIVE" | md5sum -c -

python3 - "$ARCHIVE" <<'PY'
import pathlib
import sys
import tarfile

with tarfile.open(sys.argv[1], "r|gz") as handle:
    for member in handle:
        path = pathlib.PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"unsafe archive path: {member.name}")
print("Archive path-safety check passed.")
PY

mkdir -p "$DEST"
tar --extract --gzip --file "$ARCHIVE" \
  --directory "$DEST" --no-same-owner --no-same-permissions \
  CASF-2016/coreset \
  CASF-2016/decoys_docking \
  CASF-2016/power_docking \
  CASF-2016/power_scoring \
  CASF-2016/power_ranking \
  CASF-2016/power_screening \
  CASF-2016/README.md \
  CASF-2016/benchmark_preparation.py

echo "Native ligand SDF files:"
find "$DEST/CASF-2016/coreset" -type f -name '*_ligand.sdf' | wc -l
echo "Docking-decoy MOL2 files:"
find "$DEST/CASF-2016/decoys_docking" -type f -name '*.mol2' | wc -l
echo "Extracted files:"
find "$DEST" -type f | wc -l
echo "Extracted size:"
du -sh "$DEST"
echo "Extraction complete: $DEST/CASF-2016"
