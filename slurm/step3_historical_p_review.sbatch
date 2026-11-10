#!/usr/bin/env bash
# Read-only selected historical P metadata; no rendering and no observed fits.
#SBATCH --job-name=HECTOR_HISTP
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/step3_historical_p_%j.out
#SBATCH --error=outputs/step3_historical_p_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
source env/friday_r.sh
python3 -u python/phasea/step3_historical_p_review.py \
  --repo "$PWD" --rscript "$BD_CONDA_PREFIX/bin/Rscript" \
  --out "outputs/step3_historical_p_${SLURM_JOB_ID}"
