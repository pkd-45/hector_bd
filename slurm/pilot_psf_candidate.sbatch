#!/usr/bin/env bash
# Requested PSF star fits only. Original E1-E7/G5 criteria remain pending.
#SBATCH --job-name=BD_PSF_CANDIDATE
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/pilot_psf_candidate_%j.out
#SBATCH --error=outputs/pilot_psf_candidate_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
: "${1:?Pass the validated outputs/pilot_inputs_JOB directory}"
source env/friday_r.sh
Rscript --vanilla scripts/build_pilot_psf_candidate.R \
  "$1" "$PWD/outputs/pilot_psf_candidate_${SLURM_JOB_ID}"
