# Hector/WAVES bulge-disc Track B V02 — amended Step 0 + fine-sampling review candidate

**Status (31 October 2026): SYNTHETIC-ONLY CODE / UNRUN R GATES / BLOCKED FOR REAL FITTING.** This package does not import real WAVES files, run a fit, choose a WAVES segmentation product, or modify historical directories. The authoritative specifications are `CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md` (v1.1) and `CLAUDE_BD_RESTART_03_CODE_INVENTORY.md`; Claude's V02 review approved starting synthetic-only Track B in parallel with Phase A.

## Modules

- D2 `src/segment_contract.R`: original-grid L, target labels A_g, T_g, usable Q, F_g=T_g&Q; full exclusion ledger with one-based R row/column and zero-based FITS coordinates; image/sigma/segmentation/mask hashes; fail-closed original-grid G2. **Step 0 values require evidence, while A2 permits `as_delivered_unclassified` for an as-delivered segmap. A1 adds a verified segmap-WCS position label resolver with hard stops on sky, boundary or catalogue/group contradiction.** This is array-level science verification; original FITS file/HDU provenance must be separately proven.
- D3 `src/build_inputs.R`: includes all target pixels including disconnected islands and holes, plus padding max(ceil(PSF_size/2),ceil(3*FWHM)); conserves original source coordinates and G2; supplies explicit region and native L crop. No autoclip, region redefinition, sky subtraction, segmentation modification, or likelihood fit. Sigma is the **verified provided sigma**; model-dependent shot-noise via `g_eff` is deliberately not inserted until Step 0 and Phase D decide the source units/gain. A separate sky diagnostic will be D4, not this stamp.
- D5 `src/profuse_psf_v14.R`: generated with a fingerprinted patch from the complete, byte-verified `profuse_psf_v13.R` (SHA-256 `1de310c015d027c25b90f0898cfc7011e0ab8ef3e30d889835db2d3dd569e9f0`). Removes five unnecessary MAD multipliers (one clipping, two pairs of star-sigma fallback), replaces both linear PA clip call sites with axial/circular PA clipping. Retains 21-pixel floor. Records *approximate analytic* elliptical-Moffat encircled energy within the zero-corner circle and loss for analytic/default PSFs. **This is not directly measured empirical star flux.** Empirical/apodised PSFs explicitly say EE UNMEASURED and block pilot use. The source includes `src/psf_v14_patch_source.py` for deterministic reproducibility. `src/psf_compare_v13_v14.R` compares selected star IDs and kinds on identically sourced outputs when Friday is available; no changes to real stars are claimed today.
- D6 `src/renderer_tests.R`: synthetic F1–F9 gate definitions, using real ProFit rendering and `profitLikeModel` for likelihood checks. **F2 code implemented but UNRUN:** `src/f2_finesample_test.R` passes the same SHA-256-pinned 21-pixel synthetic empirical Moffat PSF to pinned ProFit `profitSetupData(psffinesampled=FALSE)` at finesample 1, 3, 5. It renders compact n=4/Re=1.5 px and point+exponential cases using `profitLikeModel`, reports all pairwise normalized model/χ² changes and one-parameter proxy scans; it gates only 3→5 `max|Δmodel|/σ<0.1`, Δχ²<1, and proxy parameter stability. 1→3 is diagnostic for **potential** historical sensitivity only, not a measured production bias. A production adopted PSF must still be tested before pilot fitting. F1/F3/F4/F5/F6/F7/F8/F9 are unexecuted on the local machine because R is unavailable. No claim that G4 has passed.
- D7 `src/fit_runner_D7_DRAFT.R`: preflight-only design for four models S/P/D_FIXED4/BD_FREE_N, explicitly checks gates and original hashes. There is NO optimisation, fitting, MCMC or executable production runner.

## Validation, prerequisites and gate status

Local: 16/16 D1 Python synthetic tests (in separate D1 V03 ZIP), 4/4 Track B Python static guards, `bash -n` PASS, byte-verified recovered PSF v13, structural and exact-site code checks PASS. **These are not R parse or R execution tests.** Neither `Rscript` nor production FITS/RDS files are available in this sandbox.

Run on Friday in the existing pinned Conda/R environment, from a new unpacked working directory:

```bash
bash run_synthetic_tests.sh
```

This script parses all R sources first, runs R tests for D2/D3 and D5 using **synthetic matrices only**, then runs the F1–F9 renderer harness against the installed ProFit. F2 now invokes the real ProFit fine-sampling setup; if the pinned ProFit API differs or the test does not converge, the R harness fails closed and reports G4 BLOCKED. Since R is not installed here, F2 is *unexecuted*, not passed. The pinned fixture is in `tests/fixtures/f2_moffat_21px.csv`, hashed and verified at runtime. Its synthetic properties do not assert performance for real per-galaxy PSFs. No real-data paths or credentials are used.

## Evidence still required

1. Step 0 as amended by A1/A2: exact delivered WAVES segmap file/HDU/hash (subtype non-blocking), validated segmap WCS and catalogue position/label cross-check (no fallback; group members verified when present), science FITS WCS/units/sky, sigma source and physical gain.
2. Original native-array orientation versus the FITS loader: R matrix row/column is used; confirm the FITS read produces the original image's row/column mapping and verify the WCS. D3 makes no transpose decision.
3. R parse and runtime synthetic tests, then per-galaxy G2 and G3 before real input construction. F1–F9 and PSF E1–E5 must pass; for F2, native empirical PSF upsampling is performed by pinned ProFit setup, and native-grid σ/objective changes plus parameter proxy scans must pass. Later repeat with each adopted real PSF.
4. Same-source v13-v14 star IDs, selected kind and per-PSF EE; real LOOCV and PSF-spread tests for the 20-galaxy pilot. F5 analytic approximations are not a substitute for empirical measurement.
5. Frozen model details for G6: parameter bounds, priors, D_FIXED4 morphology and k, plus all model-specific constraints. Nothing is silently inherited from old R54/R61/R74 selection layers.

**Do not run any real-data fitting or the full 832 with this package.** The Phase A exporter review runs separately; it does not control the new-input segmentation definition.

## Amendment V02 / operational note

- **A1** is implemented in `bd_resolve_original_labels`. The caller must supply a function converting the catalogue RA/Dec through the *original* segmap WCS to one-based R matrix indices; the converter must be explicitly marked `verified_original_segmap_wcs=TRUE`. It stops for `L=0`, 3×3 boundaries (including adjacent background), off-frame values, inconsistent catalogue `segID`, or missing/inconsistent group membership, rather than choosing a nearby label. The 3×3 rule is deliberately conservative. Centroid offsets and magnitude differences are diagnostics only.
- **A2** accepts `step0$segmentation_product='as_delivered_unclassified'`. The actual original segmap path, HDU, digest and registered grid remain mandatory evidence; this field only records missing product-subtype documentation, never changes segmentation.
- **A3** prints `F2_RUNTIME`, `F2_METRIC` and `F2_PARAMETER_PROXY` rows. The deterministic synthetic 21×21 Moffat fixture is a **pinned test PSF**, not any fitted galaxy's observed PSF. `profitSetupData` must perform the interpolation for each fine factor, and `profitLikeModel` must return the post-convolution native-grid model. Missing API support is a failure, never a silently substituted renderer. No observed fitting or parameter optimisation is called.

**Do not run F2 on real data with this package**; the observed pilot PSF and source-noise checks are additional gates. `run_synthetic_tests.sh` is the only test entry point. Full Phase A remains independent and need not finish before the pilot, but the pilot still waits for per-object Step 0/G2 and F1–F9 execution, plus PSF gates.

**Log capture on Friday:** from the extracted `TRACK_B_V02` directory in the pinned Conda environment, run `bash run_synthetic_tests.sh 2>&1 | tee trackB_V02_synthetic.log` (use `set -o pipefail` first if running an interactive shell without it); inspect the exit status. `F2_RUNTIME` and `F2_PROFIT_BINARY_SHA256` fingerprint the actual local ProFit shared library. No source or historical data paths are modified by these tests.

## Friday order and D1 compatibility

1. **Stage 0 collector**: run the unchanged, Claude-approved D1 V04 collector from the Mac (`run_hector_bd_restart_stage0.sh friday`). The attached amendments do **not** modify that read-only collector or its historical-file safety rules.
2. **D1 V04 3-galaxy dry run**: follow that ZIP's `README_V04.md`, retaining `profound_presence_inventory.tsv`. Do not launch 13-galaxy/full-832 runs yet.
3. **Track B V02 R synthetic tests**: execute `bash run_synthetic_tests.sh 2>&1 | tee trackB_V02_synthetic.log` in this package's root inside the pinned Conda environment, and check its exit code. This calls F1–F9, including implemented F2, with synthetic arrays only.
4. Return logs, path list, array presence inventory, three-galaxy lineage tables, and this V02 package to Claude. Real 20-object fitting remains gated on file-verified Step 0 and per-object G2, successful F1–F9, and real-pilot PSF gates. No observed fitting is included here.

The attached `STEP0_DECISION_MEMO_AMENDED.md` records author-supplied A1/A2/A3 corrections that supersede the corresponding lines in the original v1.1 specification and the memo bundled with D1 V04. It does not change the original historical evidence or relabel REPORTED findings as VERIFIED.
