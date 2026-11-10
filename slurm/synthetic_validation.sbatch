#!/usr/bin/env bash
# Existing synthetic recovery and notebook/CLI parity only; no F2 or observed fit.
#SBATCH --job-name=BD_SYNTHETIC
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/synthetic_validation_%j.out
#SBATCH --error=outputs/synthetic_validation_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
BD_SOURCE="$SLURM_SUBMIT_DIR"
cd "$BD_SOURCE"
test -f .hector_bd_root
source env/friday_r.sh
export BD_SYNTHETIC_REPORT_DIR="$BD_SOURCE/outputs/synthetic_validation_${SLURM_JOB_ID}"
mkdir "$BD_SYNTHETIC_REPORT_DIR"
export TMPDIR="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/bd_synthetic_${SLURM_JOB_ID}_XXXXXX")"
trap 'rm -rf -- "$TMPDIR"' EXIT
git clone --shared --no-hardlinks --quiet "$BD_SOURCE" "$TMPDIR/source"
for BD_FILE in tests/test_synthetic_recovery.R tests/test_notebook_cli_equivalence.py tests/run_synthetic_validation.R; do
  cp "$BD_SOURCE/$BD_FILE" "$TMPDIR/source/$BD_FILE"
done
cd "$TMPDIR/source"
git rev-parse HEAD > "$BD_SYNTHETIC_REPORT_DIR/git_HEAD.txt"
git diff -- tests/test_synthetic_recovery.R tests/test_notebook_cli_equivalence.py > "$BD_SYNTHETIC_REPORT_DIR/diagnostic_changes.patch"
cp tests/run_synthetic_validation.R "$BD_SYNTHETIC_REPORT_DIR/run_synthetic_validation.R"
sha256sum tests/test_synthetic_recovery.R tests/test_notebook_cli_equivalence.py tests/run_synthetic_validation.R > "$BD_SYNTHETIC_REPORT_DIR/test_sha256.txt"
echo "JOB=$SLURM_JOB_ID NODE=$(hostname) SOURCE=$PWD R_LIBS_USER=$R_LIBS_USER"
set +e
Rscript --vanilla tests/run_synthetic_validation.R
BD_RECOVERY_STATUS=$?
python3 -u tests/test_notebook_cli_equivalence.py
BD_PARITY_STATUS=$?
set -e
printf 'RECOVERY_EXIT=%s\nPARITY_EXIT=%s\n' "$BD_RECOVERY_STATUS" "$BD_PARITY_STATUS" | tee "$BD_SYNTHETIC_REPORT_DIR/exit_status.txt"
[[ "$BD_RECOVERY_STATUS" == 0 && "$BD_PARITY_STATUS" == 0 ]]
echo SYNTHETIC_RECOVERY_AND_NOTEBOOK_CLI_PARITY_PASS
