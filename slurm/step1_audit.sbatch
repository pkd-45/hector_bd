#!/usr/bin/env bash
#SBATCH --job-name=HECTOR_SEG832
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root

: "${BD_STEP1_PYTHON:?Missing Python interpreter}"

RSCRIPT=/home/uqpdas/.conda/envs/waves_r/bin/Rscript
test -x "$RSCRIPT"
test -x "$BD_STEP1_PYTHON"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

OUT="$PWD/outputs/step1_832_${SLURM_JOB_ID}"

echo "JOB=$SLURM_JOB_ID"
echo "NODE=$(hostname)"
echo "GIT_COMMIT=$(git rev-parse HEAD)"
echo "OUTPUT=$OUT"
echo "START=$(date -Is)"

sha256sum python/phasea/step1_segmap_audit.py
"$BD_STEP1_PYTHON" -c 'import numpy, scipy, astropy'

"$BD_STEP1_PYTHON" -u python/phasea/step1_segmap_audit.py \
  --image-hdu 1 \
  --seg-hdu 1 \
  --rscript "$RSCRIPT" \
  --out "$OUT"

echo "FINISH=$(date -Is)"
cat "$OUT/summary.json"
