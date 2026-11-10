#!/usr/bin/env bash
#SBATCH --job-name=BD_PILOT_PREP
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/pilot_preparation_%j.out
#SBATCH --error=outputs/pilot_preparation_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
source env/friday_r.sh
: "${BD_STEP2_PYTHON:?Set the verified Python interpreter via sbatch --export}"
export PYTHONDONTWRITEBYTECODE=1
phase="${1:-tests}"
case "$phase" in
  tests)
    "$BD_STEP2_PYTHON" -u tests/test_freeze_pilot.py
    "$BD_STEP2_PYTHON" -u tests/test_prepare_pilot_inputs.py
    ;;
  input-tests)
    "$BD_STEP2_PYTHON" -u tests/test_prepare_pilot_inputs.py
    ;;
  candidate)
    "$BD_STEP2_PYTHON" -u scripts/freeze_pilot.py \
      --output "outputs/pilot_candidate_${SLURM_JOB_ID}.json"
    ;;
  inputs)
    "$BD_STEP2_PYTHON" -u scripts/prepare_pilot_inputs.py \
      --repo "$PWD" --out "outputs/pilot_inputs_${SLURM_JOB_ID}"
    ;;
  *) printf '%s\n' 'Expected phase tests, input-tests, candidate or inputs' >&2; exit 2 ;;
esac
