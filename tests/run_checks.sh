#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
if command -v python >/dev/null 2>&1; then
  PYTHON_BIN=python
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN=python3
else
  echo 'ERROR: Python interpreter unavailable.' >&2
  exit 2
fi
export R_LIBS_USER="$ROOT/env/restart_Rlib${R_LIBS_USER:+:$R_LIBS_USER}"
"$PYTHON_BIN" -m unittest discover -s python/phasea/tests -p 'test_*.py' -v
"$PYTHON_BIN" -m unittest discover -s tests -p 'test_*.py' -v
if ! command -v Rscript >/dev/null; then
  echo 'UNRUN: Rscript unavailable; no R parse/F1-F9/equivalence claim' >&2
  exit 2
fi
Rscript -e 'for (p in list.files("R",pattern="[.]R$",full.names=TRUE)) parse(file=p)' >/dev/null
Rscript tests/test_segment_amendment.R
Rscript tests/test_contract_synthetic.R
Rscript tests/test_sky_mask_synthetic.R
Rscript tests/test_psf_branch_contract.R
Rscript tests/test_moment_pa_convention.R
Rscript tests/test_psf_v14_synthetic.R
Rscript tests/run_renderer_tests.R
Rscript tests/test_synthetic_recovery.R
"$PYTHON_BIN" tests/test_notebook_cli_equivalence.py
