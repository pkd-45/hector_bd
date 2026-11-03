# Hector/WAVES bulge–disc Phase A D1 — review candidate V03 (R9) (31 October 2026)

**This is a read-only evidence exporter, not a galaxy-fitting pipeline.** It does not alter original WAVES segmentation, source images, FITS or historical RDS files. Historical 842-object parent and frozen 832-object R61 membership remain distinct. No R74 selection joins; R54/R61 selection is not interchanged. **No full 832 run is authorised by the present validation.** This package implements Claude's D1 V01 review R1–R3 and provisional R4–R7. R9 (segmentation attribution in the historical fitted-label frame) is fixed. R8 (exact ProFuse autoclip threshold re-evaluation) is deferred pending the actual executed historical helper/source and Friday arrays.

## Evidence currently available

- The attached reviewer test (`tests/test_lineage_causes.py`) has its original SHA-256 and is unmodified. The previous V01 produced 2/4 failing adversarial cases (island filter and autoclip). **V03 passes 4/4**, both unchanged V02-review R9 regression tests, and all 10 additional in-package tests: 16/16 total.
- V03 classifies pixel omissions against the **registered pre-tightcrop Found2Fit input domain**, never the saved tightcrop stamp rectangle. If neither the saved ProFound domain nor the registered inputs-RDS proxy establishes a boundary, missing pixels outside the stamp are **unexplained**, not called cutbox.
- The new R extractor records whether `fit$profound$image`, `segim`, `segim_orig`, `mask` exist and exports their 2D arrays if saved; absence and extraction failure are distinguished. **This still requires testing on Friday's real RDS bytes.**
- R9: missing pixels are checked against the labels actually used under the historical fitted `Data$region`, not against a majority WAVES-label mapping. If `fit$profound$segim` was not archived, available in-stamp `Data$segim` is still used; missing pixels outside the stamp without a parent segmentation remain unexplained.
- Original WAVES segmentation and `Data$segim` (and `profound$segim`, when registered) are compared as partitions, allowing one-to-one relabelling while fixing background label zero. Reported mismatch pixels are a directional modal-label discordance diagnostic, not a reconstruction of historical modified pixels.
- Exact unshifted image-copy statistics are reported separately from correlation-based registration. A uniquely matching exact ProFound segmentation subwindow can corroborate Data's offset; a conflict with verified image registration blocks the corresponding family.
- The code inspects pinned `/data/uqpdas/WAVES_v38_profit272_pilot/Rlib` explicitly, as well as current libraries. It inventories `r_1333_v*` input folders rather than silently choosing `r_1333_v25`; multiple matches block that galaxy. Non-1333 original FITS grids generate a warning, **not exclusion**.
- **No Rscript or Friday's actual RDS/FITS files exist in this sandbox**; R syntax/field compatibility, real provenance and runtime remain unverified. NumPy + synthetic adversarial cases and Bash/Python syntax are tested locally.

## Do this first: approved stage-0 collector, then three-galaxy RDS probe

On the Mac where the SSH alias `friday` works:

```bash
bash run_hector_bd_restart_stage0.sh friday
```

The script produces a new `~/Downloads/HectorBD_RESTART_STAGE0.XXXXXXXX/FRIDAY_STAGE0_EVIDENCE.txt` locally. It reads Friday evidence and historical code but makes no changes to the historical project. Transfer this ZIP to Friday and extract it in a user-owned working directory (not in historical `/data` or `/scratch` trees). In Friday's existing Conda environment:

```bash
bash src/run_phasea.sh --dry-run \
  --out "$HOME/Downloads/HectorBD_PHASEA_DRY_V03_$(date -u +%Y%m%dT%H%M%SZ)" \
  --keep-intermediates
```

A brand-new output directory is mandatory. Default `--root` is the historical FIXED4 842 parent tree and `--sample-ids` is the frozen 832 R61 roster; `--dry-run` only reads the **first three** IDs. Override a path only with the recorded actual path. When labels and segmentation/image HDUs remain undecided, the script still extracts RDS fields and produces a **three-row blocked** report plus `profound_presence_inventory.tsv`. It must not pretend to have counted original WAVES pixel omissions yet. Do not invent label IDs/HDU numbers to force a comparison.

**Send Claude and ChatGPT:** the stage-0 evidence text, `profound_presence_inventory.tsv`, `phase_a_832_rows.tsv`, `summary.json`, `runtime.json`, `r_parse_gate.txt`, the RDS `scratch` extraction metadata for the three IDs (`--keep-intermediates`), `SHA256SUMS.txt`, and any run log/stdout. The question for the first run is whether saved `fit$profound` arrays are actually present and fit the expected coordinate conventions.

## After reviewing dry-run evidence; NOT authorised yet

Once Step 0's label map, original image/segmentation product and FITS HDUs are independently established, a validated ten-pilot extension can be run (the option verifies all ten are present in the frozen 832 roster):

```bash
bash src/run_phasea.sh --pilot-13 \
 --out "$HOME/Downloads/HectorBD_PHASEA_PILOT_V03_$(date -u +%Y%m%dT%H%M%SZ)" \
 --labels "$VERIFIED_WAVES_LABEL_MAP_TSV" \
 --seg-hdu "$VERIFIED_SEG_HDU" --image-hdu "$VERIFIED_IMAGE_HDU" \
 --seg-product "$VERIFIED_WAVES_SEG_PRODUCT"
```

The name `--pilot-13` means first three + the ten frozen pilot IDs; duplicates are removed if present. All four `VERIFIED_...` variables must be supported by actual WAVES evidence. Do **not** run the full 832 yet. Any unresolved input grid, ambiguous crop coordinate or missing metadata yields a blocked family/row; it must not be promoted to a G2 pass. `Data$calcregion` is reported separately from the actual likelihood region.

## Files produced

- `profound_presence_inventory.tsv`: available even when Step 0 is not complete; per fit family and RDS field presence/dimension.
- `phase_a_832_rows.tsv`: actual row count is **3 in dry-run**, potentially 13 in pilot, and 832 only in an eventually authorised full run. Each row records requested labels, provenance, counts when defensible, source paths, PSF/sky/sigma status, support hashes, and segmentation-partition checks.
- `pixel_ledger.tsv.gz`: zero-based original FITS `(row_y0,col_x0)` and evidence-ranked omission/extra classifications. `cutbox_side_fit_crop` means outside a registered pre-tightcrop domain; with `inputs.rds` proxy only, the cause remains **inferred**, not verified as an executed cutbox call. `segim_differs_upstream` does **not** assert that an island filter specifically ran; `autoclip_NA` is observed NA evidence, not automatically proof of the executed autoclip formula.
- `overlays/`: original-seg contours, stamp boundary and missing/extra pixel positions when matching is validated.
- `run_config.json`, `runtime.json`, `summary.json`, `r_parse_gate.txt`, `SHA256SUMS.txt` and optionally `scratch/` with raw extracted RDS matrix/metadata data.

## Non-negotiable science gate

The historical support comparison is **not** G2. G2 requires the new fit's support to equal the original WAVES label set intersected with documented usable original pixels, in original FITS coordinates. Do not call `profuseFound2Fit` with defaults or perform real-data fits. Synthetic-only Track B D2/D3/D5/D6 code may proceed in parallel; no real-data inputs are authorised until Step 0 is resolved. Claude reviews the three-galaxy then pilot evidence first.
