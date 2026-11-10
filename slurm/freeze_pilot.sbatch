#!/usr/bin/env bash
# Review-only roster preparation from frozen tables; never performs a fit.
#SBATCH --job-name=BD_PILOT20
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=outputs/pilot_roster_%j.out
#SBATCH --error=outputs/pilot_roster_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through SLURM}"
cd "$SLURM_SUBMIT_DIR"
test -f .hector_bd_root
export PYTHONDONTWRITEBYTECODE=1
OUT="$PWD/outputs/pilot_roster_${SLURM_JOB_ID}"
mkdir "$OUT"
git rev-parse HEAD > "$OUT/source_commit.txt"
sha256sum scripts/freeze_pilot.py tests/test_freeze_pilot.py > "$OUT/source_sha256.txt"
python3 -u tests/test_freeze_pilot.py > "$OUT/tests.txt" 2>&1
cat "$OUT/tests.txt"
python3 -u scripts/freeze_pilot.py --output "$OUT/pilot20.json"
echo PILOT_ROSTER_PREPARED_FOR_REVIEW_NO_FIT_APPROVAL
