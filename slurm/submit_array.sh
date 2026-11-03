#!/usr/bin/env bash
set -euo pipefail
# No sbatch call in this script. It is meant to BE submitted by the user.
#SBATCH --job-name=hector_bd
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --array=1-832
# Output directory is defined at submit time: sbatch --output=... 
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
[[ -f .hector_bd_root ]] || { echo 'Missing repository marker' >&2; exit 3; }
[[ "${HECTOR_BD_FULL_APPROVED:-}" == 'YES_REVIEWED_PILOT' ]] || {
  echo 'BLOCKED: full-run manual approval environment flag missing' >&2; exit 4;
}
Rscript scripts/check_full_approval.R friday
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'SLURM_ARRAY_TASK_ID missing' >&2; exit 5; }
LIST="${HECTOR_BD_GALAXY_LIST:-}"
[[ -f "$LIST" ]] || { echo 'Declare verified 832 ID list via HECTOR_BD_GALAXY_LIST' >&2; exit 6; }
N="$(awk 'END {print NR}' "$LIST")"
[[ "$N" == '832' ]] || { echo "BLOCKED: expected 832 IDs, got $N" >&2; exit 7; }
ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")"
[[ "$ID" =~ ^W[0-9]+$ ]] || exit 8
Rscript scripts/run_one_galaxy.R --profile=friday --id="$ID" --mode=observed --run-id="${HECTOR_BD_RUN_ID:?run ID needed}"
