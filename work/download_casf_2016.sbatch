#!/bin/bash
#SBATCH --job-name=download-casf
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/home/xianyang/logs/download-casf-%j.out
#SBATCH --error=/home/xianyang/logs/download-casf-%j.out

set -euo pipefail

ROOT=/data/user_data/xianyang/casf-2016
ARCHIVE="$ROOT/casf_2016.tar.gz"
PART="$ARCHIVE.part"
EXTRACTED="$ROOT/extracted"
URL='https://zenodo.org/records/20829559/files/casf_2016.tar.gz?download=1'
EXPECTED_MD5=b433ef3da03ae94731d9bd0b67df8138

mkdir -p "$ROOT" /home/xianyang/logs

echo "Destination: $ROOT"
df -h "$ROOT"

available_bytes=$(df -PB1 "$ROOT" | awk 'NR == 2 {print $4}')
if (( available_bytes < 20000000000 )); then
  echo "ERROR: fewer than 20 GB are available; refusing to download and extract." >&2
  exit 2
fi

if [[ -f "$ARCHIVE" ]] && echo "$EXPECTED_MD5  $ARCHIVE" | md5sum -c -; then
  echo "Verified archive already exists; skipping download."
else
  curl --location --fail --show-error \
    --retry 20 --retry-delay 15 --retry-all-errors \
    --continue-at - --output "$PART" "$URL"
  mv "$PART" "$ARCHIVE"
  echo "$EXPECTED_MD5  $ARCHIVE" | md5sum -c -
fi

echo "Archive members and uncompressed bytes:"
python3 -c 'import tarfile,sys; t=tarfile.open(sys.argv[1], "r:gz"); m=t.getmembers(); print(len(m), sum(x.size for x in m))' "$ARCHIVE"

unsafe=$(tar -tzf "$ARCHIVE" | grep -E '(^/|(^|/)\.\.(/|$))' || true)
if [[ -n "$unsafe" ]]; then
  echo "ERROR: archive contains unsafe paths; refusing extraction." >&2
  echo "$unsafe" >&2
  exit 3
fi

mkdir -p "$EXTRACTED"
tar --extract --gzip --file "$ARCHIVE" \
  --directory "$EXTRACTED" --no-same-owner --no-same-permissions

echo "Extracted file count:"
find "$EXTRACTED" -type f | wc -l
echo "Pose-coordinate file count:"
find "$EXTRACTED" -type f \( \
  -iname '*.mol2' -o -iname '*.mol2.gz' -o \
  -iname '*.sdf' -o -iname '*.sdf.gz' -o \
  -iname '*.pdbqt' \
\) | wc -l
echo "Top-level extracted contents:"
find "$EXTRACTED" -mindepth 1 -maxdepth 2 -printf '%P\n' | sort | head -n 100

echo "CASF-2016 download and extraction complete."
