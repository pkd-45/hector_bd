# Hector/WAVES bulge–disc restart — V03 corrected PSF review candidate

**Status: code-delivery / synthetic verification only. No real science images have been opened, no fits have run, and no 20- or 832-object run is authorised.** This repo is designed to be dropped into an existing directory and opened as the `hector_bd/` workspace in VS Code. Every notebook and batch step calls the same functions in `R/fit_one.R` (`bd_stage`, `bd_run_one`). It deliberately fails closed on uncertain FITS and model-runtime conventions.


## V03.1 patched source (2026-11-03)

This archive **already includes** the exact reviewer-supplied final patch; do **not** apply `claude_v03_final_fixes.patch` again. A copy is preserved at `reference/reviews/claude_v03_final_fixes.patch`. Verify checksums using `sha256sum -c SHA256SUMS.txt` if your shell provides GNU `sha256sum`; on macOS use `shasum -a 256 -c SHA256SUMS.txt`. The original V03.1 manifest, verified before the local Mac setup edits, is preserved at `reference/reviews/SHA256SUMS_v03.1.txt`; the root manifest includes the setup edits.

**What you can run now:** `notebooks/00_environment_check.ipynb`, the `BD_MODE=synthetic` input/G2/renderer stages of notebook 01, and `Rscript tests/test_synthetic_recovery.R` for the full four-family synthetic fit. Notebook 01 **deliberately refuses to run its fit stage in synthetic mode**. The original Claude recovery test passed externally according to the reviewer, but this archive still contains our reconstructed fixture with the PSF correction. R execution is not available in the packaging container. The original segmap, science data, and full-832 runs remain protected.

## V03 handover: actions and verification boundary (2026-11-02)

**SOURCE:** V02 was the base. The four files referenced in the review (`claude_v02_fixes.patch`, Claude's original `test_synthetic_recovery.R`, `test_moment_pa_convention.R`, and `profit_pointsource_finesample_bug.R`) were **not attached** to the chat or available in the local source tree. Accordingly this is an independently reconstructed candidate, **not** a byte-for-byte application of Claude's patch. The original regression and bug reproducer still need to be added and executed unchanged by the reviewer.

**Implemented:** ProFit moment-PA sign correction; on-bound optimisation starts; a four-quadrant synthetic PA test; random-noise recovery fixture (original unavailable); per-start convergence, parameters, bounds and rejected-solution records; PSF-v14 analytic-Moffat parameter-model path as `modellist$psf` for all four families when verified metadata exists. An empirical-only PSF can run P **only at finesample=1** with `EMPIRICAL_P_FS1_RENDERER_LIMITATION` clearly recorded. This is a contingency, **not** a validated substitute for the analytic PSF. The observed P configuration is conservatively set to 1 pending per-object F2 acceptance.

**F2 replacement:** P uses the *analytic* ProFit PSF path at predeclared synthetic nuclear S/N=60, with 1/3/5/7 candidates compared against fine-sampling 9 and the original Δχ² < 1 and max 0.1σ thresholds. The old empirical P F2 evidence is invalidated. The compact n=4 case remains separate. P magnitude identifiability is now a gate. No F2 outcome is asserted before pinned-R execution.

**Historical D1:** `Data$finesample` and `Data$psftype` are recorded for each saved fit, including absent fields, as `historical_finesample_by_fit_json` and `historical_psftype_by_fit_json` in the 832-row table. This only reports preserved metadata, it does **not** infer historical rendering bugs from absent fields.

**Run order:** Open notebook 00 in the locally pinned `waves_r_bd` kernel and execute read-only environment checks; then run `bash tests/run_checks.sh`, including the newly reconstructed tests plus the reviewer's unchanged original recovery test. Only once those gates and Step 0 from actual Friday bytes pass should observed-mode notebook 01 be considered. Do not launch any real fit or SLURM job from this package.

**Validation available at packaging:** Python/unit/static suites, shell syntax, source hash and ZIP integrity. **No Rscript present**, so no R parse, ProFit test, F2 verdict, synthetic parameter recovery, or empirical-vs-analytic PSF reproduction can be certified here.

## V02 review handover and validation status

**No observed fit or SLURM job has been run.** V02 incorporates the actionable points from Claude's executed V01 review. Its referenced `claude_verified_fixes_R.patch` and unchanged `test_synthetic_recovery.R` were not attached to this chat, so the code fixes and recovery test have been reconstructed rather than represented as byte-identical copies. Please compare this package with Claude's original patch and run *both* recovery regressions in `waves_r`.

The shared optimizer now preflights each candidate render and fails on errors, records every L-BFGS-B convergence code and message, creates systematic PA, scale-radius and B/T starts, and adds a mandatory nested start for BD_FREE_N from the fitted D_FIXED4 solution. The moment-based PA follows the ProFit F3 convention; physical parameter recovery, model ordering and nesting are checked in the synthetic R test. The scientific fitting implementation remains an **unapproved review candidate** until that R test passes.

The segmented-image likelihood support remains the exact original WAVES `F_g`. Phase C's `sky_mask_radius_px=10` is applied **only to the independent background sampling mask**, not L, Q, F_g, or the fitting pixels. Background diagnostics record old/new sky offsets and hash the sky-only mask. Its effective-gain value must still be verified from project data.

F2 now compares 1, 3 and 5, and tests candidates against fine-sampling 9 with the same 0.1 sigma and delta chi2 < 1 constraints. The analytic-Moffat fixture check establishes PSF provenance but the high-resolution Sérsic/point-plus-disc image is a ProFit numerical reference, not an independently integrated analytic galaxy truth. `run_pilot.yaml` proposes P=1 for empirical-only fallback, other models=3, **subject to the actual F2 result and per-galaxy E1–E5**; none is preapproved. The empty `config/psf_g4_pilot.tsv` requires per-family fine-sampling evidence before observed fitting.

`tests/run_checks.sh` is the one-command gate. It executes the Phase A Python suite, project Python suite, R parsing, module tests, F1–F9, the reconstructed end-to-end synthetic recovery fit, and notebook-vs-CLI artifact identity. A missing Rscript is a nonzero **UNRUN**, never PASS. Python success alone is not full validation. The synthetic recovery test runs all four families on a known convolved bulge+disc image, checks D_FIXED4 recovery/ordering and the nested BD_FREE_N inequality, and performs no observed-data I/O.

## Mac — reuse the existing WAVES setup

The existing Mac setup uses Conda `waves_py` for Python and native R with the historical project library. It does not require a Conda environment named `waves_r`. From the directory containing the extracted `hector_bd/` folder:

```bash
cd hector_bd
conda activate waves_py
export R_LIBS_USER="$PWD/env/restart_Rlib:/Users/uqpdas/work/PhD-UQ/bt_final_v2/.Rlib"
Rscript env/install_pinned.R
Rscript -e 'if (!requireNamespace("IRkernel", quietly=TRUE)) install.packages("IRkernel", lib="env/restart_Rlib", repos="https://cloud.r-project.org")'
Rscript env/register_irkernel.R
bash tests/run_checks.sh
code .
```

Run `which Rscript` to confirm native `/usr/local/bin/Rscript` is selected. The historical library already supplies the science dependencies, `yaml` and `remotes`; eight package source SHAs differ from this restart's pins. The installer adds the exact pins in dependency order to `env/restart_Rlib`, with automatic dependency installation disabled. Any missing prerequisite is reported by R; install it into the restart library. The IRkernel command adds missing notebook tooling there as well. Launch VS Code from this shell so it inherits `R_LIBS_USER`.

The kernel ID `waves_r_bd` is retained for notebook compatibility and is registered against the R that runs `env/register_irkernel.R`; the ID does not require a Conda environment of the same name. Python dependencies for the unchanged D1 audit are `numpy`, `scipy`, `astropy`, `matplotlib` as already available in `waves_py`. Do not run the notebook on real data yet.

## Friday — up to five shell commands

```bash
cd hector_bd
conda activate waves_r
Rscript env/install_pinned.R
bash tests/run_checks.sh
Rscript env/check_environment.R
```

The `paths_friday.yaml` profile carries `/data` only as a protected root, not an invented science path. The source file index comes from the Stage 0 collector. The existing, separately reviewed D1 V04 Mac collector and read-only 3-galaxy dry run are in `python/phasea/` for traceability; their original shell safety contract is unchanged. For its invocation, use the original readme in `python/phasea/` and **not** the fitting CLI.

## VS Code / notebook

Open the workspace folder containing `.hector_bd_root`. Select kernel **R (waves_r BD pinned)**. Start at `notebooks/00_environment_check.ipynb`, then `notebooks/01_explore_one_galaxy.ipynb`. The default `BD_MODE=synthetic` prevents any fitting; `01` stops at the `fit` cell with an explicit block. Notebook stages write deterministic `*.rds` artifacts under `outputs/<run_id>/<uberID>/`; existing artifacts are read and checked against the same configuration, never silently overwritten. The CLI `scripts/run_one_galaxy.R` calls `bd_run_one`, and `bd_run_one` calls precisely the same `bd_stage` functions as the notebook cells. There is **no second fitting implementation**.

### Current gates and configuration

`config/step0.yaml` records historical settings as **REPORTED**: WAVES r-band images in microJy, AB zero point 23.9, 0.30 arcsec/pixel, upstream sky subtraction, 1333×1333 expectation, historical effective gain 687 (and old gain=1 defect). These are not VERIFIED until headers, source calibration and gain in the actual units prove them. The delivered `inputs.rds$segmap_path` is used unchanged regardless of whether its subtype is known (A2).

A1 target mapping uses the original segmap WCS to select the **exact catalogue pixel label**; `L=0` and image edge stop; a 3×3 label boundary **flags** but does not stop, and never changes A_g. A verified catalogue group mapping can add explicitly documented member labels. The exact per-object PSF gate is `config/psf_g4_pilot.tsv` (initially headers only). A synthetic F2 passing is insufficient; each actual pilot PSF needs recorded F2 and E1–E5 evidence with matching file/log SHA-256 before the fit stage unlocks. The original `T_g` and usable `F_g=T_g & Q` are compared in the original frame, preserving detached islands. No re-detection, dilation, erosion, island filter, autoclip, tightcrop, cutbox, or nearest-label fallback is allowed. The static-sky and model-dependent Poisson noise remain explicitly differentiated.

**To verify real input paths:** populate `config/verified_path_index.tsv` with the collector-confirmed paths for the pilot. HDU indices in this TSV are **zero-based FITS indices**; `io_fits.R` converts to Rfits's 1-based `ext` and forbids its extension-autoselection. The current file has headers only and therefore cannot be mistaken for verified data. Compressed `.fits.gz` input is blocked until decompression side effects are audited or a separately verified uncompressed input copy is prepared; it cannot silently decompress alongside historical files. The verified PSF RDS must be a matrix; a list-structured historical PSF requires a separately reviewed extractor, not guesswork. Detector mask policy must be `NONE_VERIFIED_FROM_SOURCE` if a separate detector mask is confirmed absent; no silent all-good assumption is permitted.

`config/run_pilot.yaml` preserves the ten named galaxies in the restart specification. The remaining ten are **not named in the source files**. `scripts/freeze_pilot.py` deterministically stratifies 10 additional IDs across size, n, and crowding **only from the frozen R61 832-object roster**, writes a JSON candidate for review and never modifies config. Do not fill with invented IDs, or mix R54/R61/R74. The actual 20 ID list and hashes must be approved and recorded before pilot execution.

### Approval sequence

1. Stage 0 collector evidence and D1 3-galaxy dry run; resolve actual files/HDUs/WCS and noise provenance.
2. D2/D3 G2 for the 20 original source maps, G3 hashes, native PSF E1–E5, and R/ProFit F1–F9 (F2 plus family-specific 5→9 accuracy included) in the pinned environment. The output is `G4_PASS` only after every actual test passes; Python static checks do not count.
3. Record all 20 approved IDs and sign `pilot_approved`, `psf_g5_approved`, `models_g6_approved` in `run_pilot.yaml`; run one observed galaxy via `01` and inspect its arrays and model card. No fit has been performed for this code delivery.
4. Claude reviews the pilot and only then authorises full-832 with `pilot_review_completed` and `full_832_approved`. `slurm/submit_array.sh` additionally requires a deliberate environment approval token; it does not itself call `sbatch`.

### Gating and scientific limitations

- **V02 R changes unexecuted in this container:** Rscript is unavailable. No claim that R source parses or that ProFit 2.8.2 model rendering converges, no F1–F9 pass, no Mac–Friday parity proof yet. `tests/test_notebook_cli_equivalence.py` is a byte-for-byte RDS parity test (synthetic stages) and must run with Rscript in `waves_r`. Existing D1 Python tests are run separately.
- **Fit backend is first-pass, explicitly gated:** L-BFGS-B deterministic multi-start for S, P, D_FIXED4 (round n=4 bulge), BD_FREE_N (11 free parameters); shared centres and model-dependent Poisson likelihood. No MCMC posterior uncertainty, BIC/DIC selection or scientifically calibrated final component priors yet. It is **not** science-approved simply because its code exists. Manual G6 includes declared bounds and optimiser behaviour.
- **Segment-radius output limitation:** the postprocess module reports intrinsic model flux within *original segment pixels* and within an **equivalent circular-radius aperture**, with stamp clipping explicitly flagged. It is not yet the fully validated Casura semi-major segment-radius estimator; a pilot acceptance gate must reject conflation. Aperture shape and model-infinite extrapolations require scientific review.
- **PSF v14:** copied from the byte-verified historical v13 patch; the larger builder depends on historical helper code not present here. Only verified matrix PSFs are consumed; this package does not rebuild astrophysical PSFs from raw stars. Existing `reference/psf_v13.R` and patch retain audit provenance.
- **Historical D1 source changed intentionally in V03:** the phasea directory preserves the approved audit code including literal protected-path security guards. Those do not select science inputs; *new* pipeline science paths are resolved solely from YAML and the verified file index.

Every output is isolated from historical data directories; no real-data fitting or full SLURM run was attempted in this delivery. Report test errors verbatim; do not mark blocked gates as passes.
