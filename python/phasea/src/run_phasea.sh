#!/usr/bin/env bash
set -euo pipefail
# Do not source historical project code or alter archived data.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ -n "${BD_PHASEA_RSCRIPT:-}" ]]; then rscript="$BD_PHASEA_RSCRIPT";
elif [[ -x /home/uqpdas/.conda/envs/waves_r/bin/Rscript ]]; then rscript=/home/uqpdas/.conda/envs/waves_r/bin/Rscript;
else rscript=$(command -v Rscript || true); fi
if [[ -z "$rscript" || ! -x "$rscript" ]]; then echo 'ERROR: Rscript missing. Set BD_PHASEA_RSCRIPT to the actual Friday Conda Rscript' >&2; exit 2; fi
py=${BD_PHASEA_PYTHON:-python3}
exec "$py" "$script_dir/export_phasea.py" --rscript "$rscript" "$@"
