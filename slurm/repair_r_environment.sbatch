#!/usr/bin/env bash
# Synthetic environment work only; never submits observed fitting.
#SBATCH --job-name=BD_R_ENV
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/r_environment_%j.out
#SBATCH --error=outputs/r_environment_%j.err

set -euo pipefail
: "${SLURM_JOB_ID:?Submit this script through sbatch}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
source env/friday_r.sh
BD_PHASE="${1:-install}"
BD_RUN="$BD_ROOT/outputs/r_environment_${SLURM_JOB_ID}"
mkdir -p "$BD_RUN"
export BD_ENV_REPORT_DIR="$BD_RUN"
export TMPDIR="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/hector_r_${SLURM_JOB_ID}_XXXXXX")"
trap 'rm -rf -- "$TMPDIR"' EXIT
echo "JOB=$SLURM_JOB_ID PHASE=$BD_PHASE NODE=$(hostname)"
echo "R=$(command -v Rscript) LIB=$R_LIBS_USER"
echo "REPORT_DIR=$BD_RUN"
sha256sum config/package_pins.yaml "$BD_CONDA_PREFIX/lib/R/etc/Makeconf"
Rscript --version

case "$BD_PHASE" in
  install)
    for BD_TOOL in "$CC" "$CXX" "$FC" "$AR" "$RANLIB"; do
      "$BD_TOOL" --version | head -n 1
    done
    Rscript --vanilla env/check_toolchain.R
    Rscript --vanilla env/install_pinned.R
    Rscript --vanilla env/verify_pinned.R
    ;;
  verify)
    Rscript --vanilla env/verify_pinned.R
    ;;
  f6)
    Rscript --vanilla tests/run_f6_test.R
    ;;
  renderers)
    # A fresh isolated F6 must pass immediately before the full suite.
    Rscript --vanilla tests/run_f6_test.R
    Rscript --vanilla tests/run_renderer_tests.R
    ;;
  *) echo "Unknown environment phase: $BD_PHASE" >&2; exit 2 ;;
esac
echo "PHASE_COMPLETED=$BD_PHASE"
