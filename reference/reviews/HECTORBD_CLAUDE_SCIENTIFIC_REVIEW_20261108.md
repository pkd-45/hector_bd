# Claude independent scientific and code audit — Hector/WAVES bulge–disc restart

**Prepared:** 2026-11-08 (Australia/Brisbane).
**Repository:** `https://github.com/pkd-45/hector_bd`
**Snapshot to review:** **TO BE FILLED FROM FRIDAY AFTER PUSH**, using exact `git rev-parse HEAD` and `git ls-remote origin refs/heads/main`.
**Last independently observed GitHub `main` before the local science push:** `7fc5ebf3dfea6096a119e87cb6608c89921fb3bc` (environment verification commit).
**Review status:** REQUEST ONLY. This document does **not** approve an observed fit, the 20-object pilot, or the full 832-object run.

> **MANDATORY PRECONDITION:** Stop if the checked-out commit differs from the newly pushed GitHub `main` or if intended scientific source changes remain uncommitted. Audit the **actual pushed source and corresponding SLURM logs**. The work reported below is not a substitute for inspecting those artifacts.

## A. What Claude must do

Perform an **independent, adversarial science and code review** of the latest Friday changes relative to `7fc5ebf3dfea6096a119e87cb6608c89921fb3bc`. Check both code and the real numerical evidence. Recompute summaries where possible, identify bugs with exact file/line/fixture, and give only the minimum safe next actions needed to establish a justified first observed-galaxy fit. Review all code changed or added, including files omitted from the agent's UI summary. Do not produce another broad historical handover.

**Science contract:** The untouched original WAVES segmentation map `L` defines target membership `T_g = (L in A_g)` for the documented label(s) `A_g`; usable-pixel indicator `Q` is determined by finite image, positive finite sigma and documented detector masks; the likelihood support is exactly `F_g = T_g & Q` in the **original pixel coordinates**. Neither new detections, island pruning, dilation, merging, hole filling, nor any per-model clipping may modify this. Independent comparison against delivered source arrays is essential. See the primary Claude plan attached to the review (`CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md`).

## B. Source priority and verification labels

Read in this order:

1. `CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md` (Claude, 2026-10-30), **the original E1–E7 and G0–G8 definitions**. It is a *plan*, not an observed-pass certificate.
2. `CLAUDE_BD_RESTART_02_FILE_LOCATIONS.md` and `reference/calibration_evidence_gap.md` (if present in the reviewed commit).
3. Actual pushed repository code (`R/`, `python/phasea/`, `scripts/`, `slurm/`, `config/`, `env/`, `tests/`, `notebooks/`) at the exact snapshot SHA.
4. Friday SLURM outputs and result manifests, including those listed below.
5. Historical read-only primary R53/R61 and FITS/RDS records when necessary. Do not modify them.

Label every statement **VERIFIED_FROM_BYTES**, **VERIFIED_FROM_CODE**, **AGENT_REPORTED**, **INFERENCE**, **UNRESOLVED**, or **NOT_RUN**. A job status `COMPLETED` alone never establishes scientific acceptance. Never elevate an agent summary to independently verified evidence.

## C. Newly reported progress to audit (not pre-approved)

### C1. Pinned runtime and renderer

The Friday agent reports:

- Exactly 27 pinned packages in `~/hector_bd/env/restart_Rlib`, manifest `outputs/r_environment_1443187/package_manifest.tsv`.
- R 4.3.3, ProFit 2.8.2, Git source SHA `2dca6ebb6b1a47ecf6ac0cb955424fc680c8bdf7`.
- `ProFit.so` SHA256 `79203e3d3055193e3009b07ec5de08f296bf0b39887a1f2c3cf3b4f3d7651b81`.
- Jobs `1443187` package/runtime PASS, `1443186` isolated F6 PASS, `1443677` full F1–F9 PASS (reported 3:24:41).
- F6 initial / post-first-render / post-100-render PSF hashes equal: `bc18ec42f8d154d3bac41bb280a740af7e6caf96a4af525923b4bc18f3b15bcd`; maximum pixel difference 0; no attribute mutation.
- F2 synthetic choices under **pinned** ProFit: 3 for compact Sérsic and 1 for analytic-Moffat point + exponential. Verify the actual family pass/accuracy/parameter criteria and runtime provenance from `outputs/r_environment_1443677.out`.

**Review questions:** Verify exact loaded library path and `RemoteSha`, all 27 package versions/source fingerprints, F6 test that compares the *caller's* actual pre/post PSF object and attributes, and whether the nine F tests cover the fitter's actual model code paths (analytic vs empirical P; sample-specific F2 gating). Do not claim global PSF validation from one synthetic PSF test.

### C2. Other synthetic tests — new job 1453593

Agent reports both `tests/test_synthetic_recovery.R` and `tests/test_notebook_cli_equivalence.py` passed in SLURM `1453593`: four fitted families converged, zero failed evaluations, eight physical-recovery assertions, BD_FREE_N objective better than D_FIXED4, byte-identical notebook/CLI stage artifacts.

**Audit:** Retrieve job stdout/stderr, exact input seed and library path, recovered parameters vs planted truth, objective comparisons, convergence codes and all physical tolerances; confirm the test exercises the current `R/fit_one.R`. Check byte parity is genuinely comparing two independently generated identical *input and G2* products, not accidentally both reading the same cached output. If the test does not cover an observed-model execution branch, explicitly say so. Treat synthetic recovery as necessary, not sufficient for observed data.

### C3. Three-galaxy sky/gain diagnostics — job 1437317

Read `outputs/step2_calibration_1437317/calibration.tsv`, `.json`, source hashes and script `python/phasea/step2_calibration.py`. Reported image BUNIT is blank; median residual sky pedestal near zero, MAD ~0.0144–0.0150 native units; coarse effective-gain estimators varied greatly: ~1162, ~784, ~406 for IDs `11690092907692`, `11690105411641`, `41691060304604`. These are **diagnostic values**, **not** verified gain or image units. The first object shows an extended positive raw ACF baseline (~0.1 at many offsets). Do not turn such correlations into BIC `N_eff` without validated treatment of low-frequency sky structure.

New `python/phasea/step2_spatial_calibration.py` reports that a **40-pixel sky exclusion** leaves inadequate spatial coverage to fit a trustworthy plane. Audit its requirements, masks, failure state and uncertainty reporting. A failed coverage criterion must be marked **UNMEASURABLE**, not silently passed, fitted anyway, or made easier until it passes. Report valid raw-noise diagnostics separately. Inspect `python/phasea/tests/test_step2_spatial_calibration.py` and the actual result of the current spatial job (if one has been run).

**Remaining primary-source blocker:** WAVES r-band tile producer/calibration metadata for the 0.30 arcsec products: flux convention (μJy? counts? rate?), zeropoint, image scaling, exposure/weight map, coadd/resampling kernel and correlation, previous sky handling, detector conversion/effective-gain definition. Neither an assumed zeropoint 23.9 nor a historical gain 687 is certified by an empty BUNIT or the rough estimator. The production documentation is attributed in the original plan to Bellstedt et al. (in preparation); an authoritative local path or public URL is **not established**. Prepare a precise short question list for the producer if unresolved. Do not manufacture the source.

### C4. Historical point-source metadata

Independently observed R61 832-row model-selection table:

`/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_R53_NATIVE_SPD_FIXED4_POLICY_EXECUTION_20260930/results/R61_NATIVE_MODEL_SELECTION_832.tsv`

Selected counts: **352 S, 68 P, 412 D_FIXED4** (=832). Native selection classes: **771 NATIVE_DECISIVE, 53 NATIVE_AMBIGUOUS, 8 NATIVE_SUPPORTED_START_SENSITIVE** (=832). These are different categorizations.

The new agent reports all **68 selected P** map uniquely to recorded R53 run and start IDs; the field extraction was queued, **not confirmed complete**. Audit `python/phasea/src/export_historical_p_metadata.R` and actual R53 file resolution, and inspect extracted `Data$finesample`, `Data$psftype`, saved PSF source, `ProFit` version/ABI where available, explicit missing-field flags, exact selected-run provenance and selection leakage. **Never infer a value from an absent field, the present runtime, or a fallback default.** Report counts by actual value/absent/error with source file hashes.

### C5. Original 832-input and WCS audits

Read-only SLURM Step-1 job `1427444` reported 832/832 catalogue-label agreement, 822 exact historical/delivered segmentation comparisons and 10 with a small historical target-only island removed (26 target pixels in aggregate). 24 objects have science/segmap header WCS offsets >0.05 arcsec; separate block-matching evidence historically gave 9 `MISREGISTERED_CANDIDATE`, 11 `AMBIGUOUS`, 4 `AMBIGUOUS_INSUFFICIENT`. **Do not recast candidates as confirmed WCS errors** without independent image/pixel evidence; do not shift labels. Four flagged IDs are `41680096212907`, `12690079611187`, `12690076908124`, `43690073113239`; exclude from a clean observed pilot pending adjudication, without altering the frozen 832 roster. Distinguish this preliminary 832 audit from other 842-parent and R74 sample counts. Frozen roster original SHA256 reported `874a350440b0361371db15b4f4e363bfc79eb924a31ff0ee5ae478ca41d275f3`.

### C6. Pilot freeze and control-flow changes

Review the newly edited `scripts/freeze_pilot.py`, new/edited `slurm/freeze_pilot.sbatch`, `slurm/pilot_preparation.sbatch`, all selection fixtures and pilot output. The old helper assumed `R50`, physical Sérsic `n` and `crowding`, absent from actual R61 table. Agent reports new strata use recorded catalogue brightness/size, resolution quartiles and R61 model classes. Audit joins, missingness, duplicate handling, deterministic seed, no outcome leakage, frozen cohort integrity, objective pre-fit selection, and appropriate representation of compact, PSF-dominated, edge-on, and ambiguous cases. Verify policy regarding historically named `W12690079611187`, itself WCS-flagged. **Do not quietly reselect a scientifically motivated known case or edit the 832 population**. Publish both a documented exception decision and a reproducible frozen list for human review.

## D. Authoritative PSF and acceptance gates

The text below is grounded in `CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md` §7 and §14; quote that file as the source of record and inspect the exact original before accepting.

| Item | Original plan's requirement; evidence expected |
|---|---|
| **E-1** | Fix double application of MAD scale (`1.4826`) in star clipping and σ fallback; rerun selection and identify changed stars. |
| **E-2** | Position-angle clipping using doubled-angle circular statistics (axial PA). |
| **E-3** | Read `d_use`, FWHM and β from saved PSFs, calculate encircled energy of adopted (including zeroed-corner) PSF versus analytic Moffat; retain ≥99% if possible or explicitly flag loss; no silent renormalization hiding missing flux. |
| **E-4** | Record PSF kind (analytic, empirical or default) per object, and use the *same* PSF kernel for fitting, scoring and display. |
| **E-5** | Verify ProFit ownership / no caller PSF mutation with unchanged hashes over repeated evaluations. |
| **E-6** | Leave-one-star-out residuals, reasonable-PSF spread, and propagation of that spread into fitted `n`, `Re`, `B/T` for the pilot. |
| **E-7** | Spatial PSF verification: quadrant-nearest stars and FWHM variation over the ~400 arcsec field. |

Gates from the original plan (all must retain their blocking nature):

- **G0**: Step 0 S0.1–S0.3 provenance recorded from evidence.
- **G1**: Complete 832-object Phase A pixel lineage; every missing pixel classified.
- **G2**: `F_g == T_g & Q` on the original coordinate grid, with exclusion ledger.
- **G3**: Input and model-relevant source hashes unchanged across stages.
- **G4**: F1–F9 all passed under pinned runtime (reported passed; verify).
- **G5**: PSF E1–E5 resolved, with E6 reported for pilot.
- **G6**: Predeclare models, bounds, parameter count k and prior policy.
- **G7**: Pilot cards reviewed; no model uses pixels outside original `F_g` (before full campaign).
- **G8**: Model selection calibrated and validated on separate subsets before science claims.

Do not confuse F2 family-level synthetic convergence with per-object E1–E5 PSF acceptance.

## E. Independent adversarial review questions

1. **Provenance and exactness:** Are all intended source changes online, reproducible and included in the reviewed SHA? Compare working-tree status, `origin/main`, committed diff, manifests and SLURM script hashes. Flag any local-only files, stale run evidence, source/log mismatch or missing input indexes.
2. **New code bugs:** Inspect each new or changed Python/R/sbatch file. Check dimensions, FITS HDUs and array orientation, WCS/grid conventions, segment membership, NaN/mask handling, nonfinite values, implicit fallback, multiprocessing/SLURM environment, unsafe output overwrites and unintended data writes. Give tests and minimal patches for genuine defects.
3. **Sky/noise science:** Does the estimator use independent sky pixels? Is the nominal sky-mask radius defensible? Are flat gradients and background residuals measurable or honestly unmeasurable? Are MAD, normalized residuals and ACF correctly computed? Does correlated noise affect likelihood and information criterion; what is validated vs only explored? What authoritative source must settle gain and zeropoint?
4. **PSF:** Are E1–E5 actual observed-data evidence rather than synthetic-only proofs? Is F6 fixing a real mutation of caller data? Are analytic and empirical branches consistent and traceable to saved PSF provenance? Is E6/E7 evidence planned correctly?
5. **Historical P:** For all 68 selected P records, verify exact R53 path/run/start joins and reports of present, absent and unparsable finesample/psftype. No default-value filling or model-selection reclassification.
6. **Pilot design:** Does the revised scheme select 20 unique eligible IDs without post-hoc cherry-picking? Are known WCS exceptions handled openly? Does it preserve a useful scientific control sample? Is the script's output immutable and hashable after review?
7. **Synthetic tests:** Are recovery assertions meaningful and independent? Did byte parity test include only synthetic stages, not observed modes? Is nesting BD_FREE_N vs D_FIXED4 valid under the same support and PSF? Any objective/likelihood sign errors or hidden failure plateaus?
8. **Approval barriers:** Demonstrate that a real fit remains blocked without Step 0 evidence, per-object PSF index, configured approvals and correct input hashes; no full-run bypass. Negative tests should explicitly try to violate these gates.
9. **Alternative scientific interpretations:** Could the 3-gal gain spread reflect flux units, varying exposure weights, galaxy contamination, correlated warping or model-dependent noise assumptions rather than literal detector gain variation? Identify discriminating measurements, not guesses.
10. **Efficiency:** Avoid resubmitting an already completed 3.5-hour F2 job or recalculating 832 source arrays unless specifically justified by a confirmed defect. Read the existing outputs first.

## F. Expected Claude deliverables

Return ONE cohesive Markdown report, with:

1. **Verdict table** for: Git publish integrity, pinned R runtime, F1–F9, synthetic recovery, notebook parity, 3-gal calibration, spatial coverage, historical P provenance, pilot freeze, Step 0, and readiness for *one observed fit*; use **PASS / PASS_WITH_LIMITS / FAIL / BLOCKED / NOT_RUN** with source path/line and job proof.
2. **Issue register**, sorted P0 (can corrupt science, invalid data lineage or bypass approvals), P1 (blocks defensible inference), P2 (quality/efficiency); each issue: file:line, reproducible evidence, minimal fix, regression test.
3. **Numerical rechecks** for all asserted counts (27 packages, 9 renderer tests, 8 recovery assertions, 68 P, 832 sample, 24 WCS, 10 islands). Explicitly separate independently recomputed numbers from agent reports.
4. **Evidence ledger**: source repository SHA, log file paths and SHA256, actual native inputs consulted, the test environment and exact commands.
5. **Go/no-go decision** for first observed `11690092907692`, for 20-object pilot, and for 832 full run **separately**. Do not infer approval from synthetic PASS.
6. **Shortest actionable next sequence**: no more than 3 scientific steps plus any P0 repairs. Supply runnable changes where possible; compute exclusively via SLURM and no silent job throttles.
7. Explicit **questions for WAVES tile producers** if authoritative producer metadata remains unavailable, constrained to exactly what's necessary to certify units, zeropoint, sky, exposure-weight/gain and covariance.

## G. Hard safety limits

- **Never fit real galaxies merely because F1–F9/synthetic recovery pass.** First observed galaxy and pilot require explicit prior approval after Step 0/G5/G6 evidence.
- Do not launch a full 832 fit, rewrite segmentation, re-cut/resample images, relax coverage limits to obtain a PASS, shift WCS or change the cohort silently.
- All heavy work on Friday must run as SLURM jobs, not on the login node.
- Historical `/data` and `/scratch` source inputs and the original model products are read-only.
- Preserve package pins and the verified ProFit ownership fix; any proposed code change must have a reason and a regression test.
- Do not invent missing files or science output values. If the actual Friday artifacts are inaccessible, write **NOT INSPECTED** and state precisely what to obtain.

## H. Evidence manifest — complete immediately before requesting Claude review

The person/agent publishing this review must replace the following placeholders from Friday **after `git push`**, or include a separate generated provenance manifest in the review bundle:

```text
PUSHED_MAIN_SHA=REQUIRED
FRIDAY_HEAD_SHA=REQUIRED
LOCAL_WORKTREE_CLEAN=REQUIRED_TRUE
REMOTE_HEAD_EQUALS_LOCAL=REQUIRED_TRUE
REVIEW_BASE_SHA=7fc5ebf3dfea6096a119e87cb6608c89921fb3bc
COMMITTED_DIFF=git diff --stat REVIEW_BASE_SHA..HEAD
SLURM_1453593_STDOUT_PATH=REQUIRED
SLURM_1453593_STDERR_PATH=REQUIRED
P_METADATA_68_RESULT_PATH=IF_COMPLETE_OTHERWISE_NOT_RUN
PILOT_FROZEN_RESULT_PATH=IF_COMPLETE_OTHERWISE_NOT_RUN
SPATIAL_CALIBRATION_RESULT_PATH=IF_COMPLETE_OTHERWISE_NOT_RUN
TILE_PRODUCER_METADATA=UNRESOLVED_UNLESS_PRIMARY_SOURCE_FOUND
FIRST_OBSERVED_FIT=NOT_AUTHORIZED
FULL_832_RUN=NOT_AUTHORIZED
```

**Review instruction:** Independent scientific audit, not affirmation. If any source is missing, keep the appropriate gate blocked. End with a concrete accept/reject/fix list and an exact next command only when scientifically safe.
