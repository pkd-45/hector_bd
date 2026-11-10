#!/usr/bin/env bash
#SBATCH --job-name=HECTOR_SPATIAL3
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/step2_spatial_%j.out
#SBATCH --error=outputs/step2_spatial_%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
: "${BD_STEP2_PYTHON:?Set the verified Python interpreter via sbatch --export}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="${TMPDIR:-/tmp}/hector_spatial_mpl_${SLURM_JOB_ID}"
phase="${1:-tests}"
case "$phase" in
  tests)
    "$BD_STEP2_PYTHON" -u python/phasea/tests/test_step2_spatial_calibration.py
    printf '%s\n' SPATIAL_CALIBRATION_SYNTHETIC_TESTS_PASS
    ;;
  measure)
    "$BD_STEP2_PYTHON" -u python/phasea/step2_spatial_calibration.py \
      --repo "$PWD" \
      --audit outputs/step1_catalogue_1427444/step1_segmap_audit.tsv \
      --out "outputs/step2_spatial_calibration_${SLURM_JOB_ID}"
    ;;
  *) printf '%s\n' 'Expected phase tests or measure' >&2; exit 2 ;;
esac
