# D1 V04 delta — original-label extra-pixel counts + Step 0 byte checks

The existing D1 V03 read-only and fail-closed contracts remain in force. This patch builds on it without authorising any real-data fit.

## Changes
- Per family and selected galaxy: `N_extra_original_background` (original WAVES L=0) and `N_extra_original_neighbour` (L nonzero, not in target labels). Both are counted in `N_extra` and each extra pixel is labelled in `pixel_ledger.tsv.gz` with its original label/provenance. Background extra is `DILATION_CONSISTENT_NOT_CAUSAL_PROOF`; neighbour extras indicate observed encroachment, not why it occurred.
- `summary.json` adds `per_extra_pixel_class_selected_only` counts. The missing-pixel classifications and R9 lineage logic are untouched.
- Stage 0 collector retains `paths.tsv` to the new output directory and prints machine-friendly reported-vs-byte-verified image/segmentation checks and `profoundGainEst` signature. No gain estimation is executed automatically.
- `STEP0_DECISION_MEMO.md` states all old settings as REPORTED, exact validation sources, and the unresolved distinction between catalogue mapping and WCS centre pick.

## Execution (on the Mac with SSH alias `friday`)
1. `bash run_hector_bd_restart_stage0.sh friday` (no historical writes; produces `~/Downloads/HectorBD_RESTART_STAGE0.*/FRIDAY_STAGE0_EVIDENCE.txt` plus `paths.tsv`).
2. Transfer this ZIP to Friday, unpack in a fresh user-owned nonhistorical directory. In its existing pinned Conda/R environment, run:
   `bash src/run_phasea.sh --dry-run --keep-intermediates --out "$HOME/Downloads/HectorBD_PHASEA_DRY_V04_$(date -u +%Y%m%dT%H%M%SZ)"`
3. Separately, from the Track B V01 package, run `bash run_synthetic_tests.sh` in the pinned R/ProFit environment; **F2 may remain explicitly BLOCKED**. Do not call it passed if blocked.
4. Upload to this chat and Claude the collector's log/`paths.tsv`, three-galaxy output directory contents including `profound_presence_inventory.tsv`, and Track B R-test output.

Do **not** run `--pilot-13`/the full 832, build real fit inputs, or fit any galaxy on the strength of this source-only test. Check the evidence/gates first. The original D1 `README.md` retains detailed arguments, run outputs, test requirements and safe-output restrictions.

## Local validation of this patch
- Python: 18 unit tests, including unchanged adversarial tests, must pass.
- Bash parse: `bash -n run_hector_bd_restart_stage0.sh && bash -n src/run_phasea.sh`.
- Python syntax: `python -m py_compile src/export_phasea.py` and embedded collector Python parsed.
- R parsing/Friday RDS and actual FITS evidence: **not available in this environment; must run on Friday**.
