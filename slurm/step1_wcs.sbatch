#!/usr/bin/env bash
#SBATCH --job-name=HECTOR_WCS24
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root

: "${BD_STEP1_PYTHON:?Python interpreter required}"
test -x "$BD_STEP1_PYTHON"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

AUDIT=outputs/step1_832_1427058/step1_segmap_audit.tsv
OUT=outputs/step1_832_1427058/wcs_registration.tsv

test -f "$AUDIT"
test -f python/phasea/wcs_registration_blockxcorr.py
test ! -e "$OUT"

echo "NODE=$(hostname)"
echo "JOB=$SLURM_JOB_ID"
echo "START=$(date -Is)"

"$BD_STEP1_PYTHON" -u \
  python/phasea/wcs_registration_blockxcorr.py \
  --audit-tsv "$AUDIT" \
  --image-hdu 1 \
  --seg-hdu 1 \
  --min-wcs-arcsec 0.05 \
  --search 15 \
  --out "$OUT"

echo "FINISH=$(date -Is)"
cat outputs/step1_832_1427058/wcs_registration.summary.json
