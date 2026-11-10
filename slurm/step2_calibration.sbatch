#!/usr/bin/env bash
#SBATCH --job-name=HECTOR_CAL3
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
: "${BD_STEP2_PYTHON:?Set Python interpreter via sbatch --export}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$BD_STEP2_PYTHON" -u python/phasea/step2_calibration.py \
  --repo "$PWD" \
  --audit outputs/step1_catalogue_1427444/step1_segmap_audit.tsv \
  --out "outputs/step2_calibration_${SLURM_JOB_ID}"
