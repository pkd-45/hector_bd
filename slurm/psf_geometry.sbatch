#!/usr/bin/env bash
#SBATCH --job-name=BD_PSF_GEOMETRY
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/psf_geometry_%j.out
#SBATCH --error=outputs/psf_geometry_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
source env/friday_r.sh
export BD_PSF_GEOMETRY_REPORT_DIR="$PWD/outputs/psf_geometry_${SLURM_JOB_ID}"
mkdir "$BD_PSF_GEOMETRY_REPORT_DIR"
Rscript --vanilla tests/test_psf_centering.R
