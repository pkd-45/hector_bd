# Changelog

## Local Mac setup correction — 2026-11-03

- Reuse `waves_py`, native R and the existing WAVES `.Rlib`; keep exact restart pins in `env/restart_Rlib`.
- Handle the atomic `NA` returned for an absent package in the pinned installer, and install pins in dependency order without fetching unpinned remote dependencies.
- Force installation of packages missing from the restart library even when `remotes` finds the same SHA in a historical library; existing matching restart packages still skip.
- Add the restart library to the test runner and select it before loading modules in the environment check; point Mac VS Code R settings at native R.
- The two V03.1 fitting fixes remain present. These setup edits do not repeat or recertify the reviewer's reported R fitting runs.

## V01 — 2026-10-31
- Converted the reviewed Track B V02 arrays, PSF v14 and F1–F9 into a standalone `hector_bd/` workspace.
- A1 correction: 3×3 label boundary flagged, exact WCS centre label used; L=0/image edge stop.
- Added unified `bd_stage()` and `bd_run_one()` notebook/CLI path, fail-closed resume and immutable SHA-256 input hashes.
- Added explicit Mac/Friday configuration profiles, unverified Step 0 evidence ledger, frozen seven package Git SHAs, and separate restart Rlib installer.
- Added observed FITS/WCS index loader with explicit zero-based HDU mapping, synthetic fixed-galaxy input, model family definitions, prototype optimiser, cards and aperture diagnostics.
- Preserved original D1 V04 Python exporter unchanged under `python/phasea` and source ProFit historical PSF diff under `reference/`.
- Added IRkernel notebooks, SLURM approval gates, R/Python tests and synthetic notebook/CLI byte equivalence test (not executed without Rscript).
- Left the ten as-yet-unselected stratified pilot IDs unresolved rather than invent them; all real-data approvals default false.

## V02 review candidate — 2026-11-02

- Reconstructed the independently executed Claude fixes from the written review. The original `claude_verified_fixes_R.patch` and unchanged `test_synthetic_recovery.R` were **not attached** in this conversation; this version does not claim to have applied or included those exact files.
- Fixed `%in%` dimension loss in WCS label-centroid computation, nested `tofit`/`tolog` flags, ProFit check rendering, F8/F9 call signature, and `profitRemakeModellist` spelling.
- F2 avoids ProFit's no-free-parameter `inheritIDs` failure; independently checks the native fixture against an analytic Moffat, evaluates high-finesample 9 reference, and applies the same delta-chi-square and noise-unit thresholds to candidate family fine-sampling factors. F2/PSF approval remains BLOCKED until verified on pinned R.
- Replaced weak Gaussian perturbations with second-moment geometry and declared radius/PA/light-split start variants. BD_FREE_N is also seeded by converged D_FIXED4 with exact nested parameters. All renders preflight, objective failures stop immediately, convergence and messages persist.
- Added a sky-only dilated exclusion mask for background diagnostics; original source-label segmentation is never modified. Expanded per-object F2 manifests to record the accepted sampling factor for all four families.
- Pinned magicaxis, celestial, ParmOff and cmaeshpc to exact upstream source SHAs, using ICRAR/ProFit for ProFit.
- Added an independent reconstructed end-to-end synthetic recovery regression, a sky-mask regression, and static guard checks to `tests/run_checks.sh`. **Rscript and pinned ProFit are unavailable in the packaging container, so the R suite and its recovery assertion remain unexecuted.** No observed fits or SLURM jobs launched.

## V03 independent review candidate — 2026-11-02

- Fixed 180°-mirror ProFit moment angle with negative `atan2` sign; four-quadrant synthetic R regression added (not executed here).
- Changed structured-start bounds clipping from inward 1e-6 to exact bounds to preserve D_FIXED4→BD_FREE_N embedding.
- Recorded every converged and uncertified start, including optimiser code/message/parameters and bound hits; emitted per-family diagnostic RDS even when all starts fail.
- Introduced explicit analytic Moffat PSF parameter-model path (`modellist$psf`, `psfdim >=81`, no image PSF) for all four model families where verified v14 metadata exist. Empirical-only P is blocked above fine-sampling 1 and flagged.
- Replaced P F2 bitmap-PSF comparison with analytic model at fixed synthetic peak S/N 60, added fs7 candidate vs fs9 reference, and retained strict Δχ² < 1, 0.1 sigma and parameter tests. No results claimed until pinned-R tests.
- Reconstructed random-noise synthetic recovery fixture (≥500 optimiser iterations); original Claude fixture/patch was not attached, so must be run separately unchanged.
- Added historical `Data$finesample` and `Data$psftype` extraction and distinct per-fit JSON columns to D1 832-row exporter.
- Real-data approvals remain false, ten stratified pilot IDs remain unfilled, all historical paths remain read-only.

## V03.1 source integration — 2026-11-03

- Applied reviewer-supplied `reference/reviews/claude_v03_final_fixes.patch` **exactly** to the V03 source using `patch -p1` with no conflicts. Only `R/fit_one.R` and `tests/test_synthetic_recovery.R` changed due to the patch.
- `bd_fit_four_models` now passes the normalised, fingerprinted `built$input$psf` into the PSF selector. The independent synthetic recovery fixture supplies its PSF explicitly.
- The reviewer reports pinned-R execution: original PA test, F1/F3–F9, F2 analytic-Moffat and original recovery, plus notebook/CLI byte identity passed. Those tests have **not been rerun** in this packaging environment (no Rscript).
- **Operational distinction:** Notebook 01 in `BD_MODE=synthetic` intentionally blocks at `bd_stage(ctx, "fit")`; it checks inputs/G2/renderer but does not run synthetic optimisation. Synthetic four-family recovery is exercised by `Rscript tests/test_synthetic_recovery.R` instead. Do not misrepresent this as an end-to-end fit from the notebook. No observed fit is authorised without file-evidenced Step 0 and explicit approvals.
