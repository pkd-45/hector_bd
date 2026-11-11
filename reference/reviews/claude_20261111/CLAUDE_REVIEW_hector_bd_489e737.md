# Independent review — `pkd-45/hector_bd` at `489e73742b9cfe0e303c4c4162cf48bf5e74526f`

**Reviewer:** Claude. **Date:** 11 November 2026.
**Environment used for recomputation:** Ubuntu 24.04 sandbox, R 4.3.3, Python 3.12; ProFit 2.8.2 built from `2dca6eb` and the other ten pinned packages built from their pinned SHAs (CRAN dependencies from Ubuntu, not the 16 CRAN locks). Observed-data arrays were **not** available to me; every observed number below is recomputed from the Friday output files in the bundle, never from an agent summary.

Labels: **VERIFIED_FROM_BYTES** (I recomputed it from files), **VERIFIED_FROM_CODE** (read in the committed source), **AGENT_REPORTED** (only the bundle says so), **INFERENCE**, **UNRESOLVED**, **NOT_RUN**.

---

## 1. Verdict table

| Item | Verdict | Evidence |
|---|---|---|
| Git publish integrity | **PASS** | VERIFIED_FROM_BYTES: fresh clone at `489e737`, `SHA256SUMS.txt` passes; `diff -r` of `PUSHED_SOURCE/` against the clone: no differences. `PUSH_PROVENANCE.md`: tree clean, local == origin. |
| Pinned R runtime (27 packages, ProFit `2dca6eb`) | **PASS** | VERIFIED_FROM_BYTES: `r_environment_1443187/package_manifest.tsv` has 27 rows; ProFit row = `2.8.2`, `ICRAR/ProFit`, `2dca6eb…`; `ProFit.so` SHA `79203e3d…` appears identically in jobs 1443677, 1453593 and 1443187. |
| F1–F9 (job 1443677) | **PASS** | VERIFIED_FROM_BYTES: nine `Fn PASS` lines; F2 `FAMILY_PASS TRUE` for both cases (compact → finesample 3; analytic-Moffat point source → 1); F6 PSF hash unchanged across 100 renders. Matches my own earlier runs of the same code. |
| Synthetic recovery (job 1453593) | **PASS_WITH_LIMITS** | VERIFIED_FROM_BYTES: four families converged (code 0), 0 failed evaluations, BD_FREE_N −1249.42 ≤ D_FIXED4 −1249.38, 8/8 physical checks. **Limit:** the tolerances are loose (bulge mag 0.65, bulge Re 1.5 px on 3.5 px, PA 15°); actual errors 0.18 mag and 15 % in bulge Re would fail my original test. Test file hashes match the commit. |
| Notebook/CLI parity | **PASS** | VERIFIED_FROM_BYTES: `inputs.rds` and `g2.rds` byte-identical (SHA `584f38…`, `700973…`), synthetic mode only, as intended. |
| PSF geometry (job 1453626, 15/21 fail) | **FAIL → diagnosed, patch supplied** | VERIFIED_FROM_CODE + reproduced: `R/psf_v14.R:152` places the Moffat at `(dim+1)/2`; ProFit's grid centre is `dim/2`. Offset → 0.5 px exactly as dim grows (0.493/0.498/0.500). Same line in historical `reference/psf_v13.R:152`. |
| 3-galaxy calibration (job 1437317) | **PASS_WITH_LIMITS** (diagnostic only) | VERIFIED_FROM_BYTES: pedestal −0.0002 to −0.0005, σ_MAD 0.0144–0.0150, blank-sky z: std 1.007–1.016, skew 0.02–0.03, |z|>3 fraction 0.36–0.43 %. Gain estimates 1162 / 784 / 406 are ProFound `profoundGainEst` lower bounds, not a measurement. BUNIT and every calibration keyword absent. |
| Spatial sky / N_eff (job 1453640) | **UNMEASURABLE for the plane; short-range ACF measured** | VERIFIED_FROM_BYTES: two of three fields share ρ(1)≈−0.075, ρ(2)≈+0.05, plateau ≈0.003 (Σρ over lags ≤5 = 1.6–2.2, N_eff ≈ N/2); `11690092907692` adds a 0.05–0.07 plateau to 64 px that plane detrending does not remove; 40 px mask → “too few well-distributed cells” (job 1453615). |
| Historical P provenance (job 1453616) | **PASS_WITH_LIMITS** | VERIFIED_FROM_BYTES: 68 unique records; all 68 `Data$finesample = 1`, `Data$psftype = 'empirical'`, `saved_ProFit_version = 2.7.2` in both layers; `Data$psffinesampled` absent in all 68 (expected for 2.7.2). **Limit:** the executed helper (SHA `7331113…`) ≠ committed helper (`5297a3c…`); `summary.json` missing. |
| Pilot freeze (job 1453629) | **PASS** | VERIFIED_FROM_BYTES: 20 unique IDs; none of the 24 WCS IDs; 33 exclusions = 24 WCS + 10 pruned − 1 overlap; `W12690079611187` → `W11690092907692` replacement recorded; coverage 7 D / 6 P / 7 S, all brightness×size strata; `config_modified: false`. |
| Rfits input copy (job 1453635) | **FAIL → diagnosed, patch supplied** | VERIFIED_FROM_BYTES: all 10 samples agree to ≤3.2×10⁻¹⁵ relative; R wrote 15 significant digits, Python compared with `!=`. Full-array hash identical. Orientation **is** established. |
| Step 0 (`config/step0.yaml`) | **BLOCKED** | VERIFIED_FROM_CODE: eight fields `REPORTED`, gain `REPORTED`; `bd_step0_gate()` (`R/fit_one.R:25–36`) correctly refuses. |
| Readiness for one observed fit | **BLOCKED** | See §6. |

---

## 2. Numerical rechecks (independently recomputed)

| Claim | Recomputed | Source |
|---|---|---|
| 27 pinned packages | 27 | `package_manifest.tsv` (28 lines incl. header) |
| 9 renderer tests pass | 9 (F6 printed twice: isolated run + suite) | `r_environment_1443677.out` |
| 8 recovery assertions | 8, all pass | `synthetic_recovery_report.json → physical_recovery` |
| 68 selected P | 68 records, 68 unique IDs | `historical_p_metadata.jsonl` |
| 832 roster; 352 S / 68 P / 412 D_FIXED4 | 832 rows; counts in `pilot20.json → native_model_counts` and `binding_evidence.json` | both files |
| 24 WCS header exceptions | 24 rows | `wcs_registration.tsv` |
| 10 island-pruned | 10 `TARGET_REMOVED_ONLY` | `step1_segmap_audit.tsv` (1427444) |
| 15/21 PSF cases fail | 15 fail centring, 15 fail symmetry; all 6 `direct_physical_centre_*` pass | `centering_metrics.tsv` |
| 13 kernel differences | 13 accepted-vs-source `DIFFERENT_NUMERIC_KERNEL_VALUES`, 55 exact; separately 60/66 accepted-vs-cache differ, 6 exact, 2 no cache | `historical_p_metadata.jsonl` |
| 7 missing items | 2 harmless glob misses + 5 missing outputs of two failed jobs | `MANIFEST_SHA256.json → missing_or_omitted` |

---

## 3. Issue register

**P0 — none.** No result in the bundle is scientifically wrong, no historical input was written, no gate can be bypassed.

### P1-1 PSF rendered half a pixel off centre — `R/psf_v14.R:152`, also `:518`, `:783`
`cx <- (dim[1] + 1) / 2` uses the 1-based index of the central pixel; ProFit's coordinate of that pixel's centre is `dim/2`. Every v14 (and v13) analytic PSF image is therefore displaced by +0.5 px in x and y. **Nature:** coordinate-convention defect in the code — not physical, not numerical precision, not a test artefact (the test's own `direct()` render at `d/2` is centred to 10⁻¹⁶).
**Reproduced:** my `/tmp/psf_centre_check.R` on the committed code gives Δ = +0.497 (21 px), +0.499 (31), +0.500 (61); with the patch Δ = 0.0000 and asymmetry 2×10⁻¹⁶ at all sizes.
**Science impact (INFERENCE):** the convolved model shifts rigidly by 0.5 px and the free centre absorbs it, so structural parameters are unbiased; fitted centres are offset by −0.5 px; kernel truncation becomes slightly asymmetric; and in v2 the P family (analytic `modellist$psf`, correctly centred by ProFit) would be inconsistent by 0.5 px with the image-PSF families. Historical v1 fits carried the same offset in the centres only.
**Fix:** `claude_psf_centre_fix.patch` (3 sites: render centre; star-stamp centre used for the initial guess and the 1 px acceptance window, which was asymmetric by 0.5 px). **Regression:** the existing `tests/test_psf_centering.R` must pass (expected after patch: all 21 cases centred); `test_psf_v14_synthetic.R` and `test_psf_branch_contract.R` still pass with the patch (run here).

### P1-2 Rfits sample comparison uses exact float equality on 15-digit JSON — `scripts/prepare_pilot_inputs.py:163–166`
`verify_pilot_fits.R:52` writes `jsonlite::toJSON(..., digits=NA)` = 15 significant digits; Python compares `expected['value'] != actual['value']`. All ten samples agree to ≤3.2×10⁻¹⁵; the comparison fails anyway and the job aborts before writing `originals_after.json`, `candidate_path_index.tsv`, `candidate_step0.json` and `SHA256SUMS.txt`.
**Orientation question answered:** the full-array hash (R column-major `[x,y]` bytes = Astropy C-order `[y,x]` bytes) **and** the ten point samples, including the off-diagonal pair (2,5)/(5,2) with different values, both establish `Rfits[x,y] == Astropy[y-1,x-1]`. A transposed read would have failed at those two points.
**Fix:** `claude_rfits_sample_tolerance.patch` (`math.isclose(rel_tol=1e-12, abs_tol=1e-15)`); `test_rfits_sample_tolerance.py` accepts the real 15-digit case and rejects a transposed read.

### P1-3 Executed historical-P helper is not the committed helper — provenance
`binding_evidence.json → helper_sha256 = 7331113…`; committed `python/phasea/src/export_historical_p_metadata.R` = `5297a3c…` (parses cleanly, 127 lines, ends at the completion message). The executed copy had a broken trailing expression (`unexpected ',' in "lags=I(unique(flags)),"`), reached only after the loop because `Rscript` evaluates top-level expressions sequentially. The 68 records are therefore valid, but **not reproducible from the committed source**, and the driver's `SystemExit` is why `summary.json` is absent.
**Fix:** none to code; rerun `slurm/step3_historical_p_review.sbatch` from the committed tree (read-only, ~10 min) so the helper hash in the evidence equals the commit and `summary.json` exists.

### P2-1 `compare_psf` cannot distinguish renormalisation from a different kernel — `export_historical_p_metadata.R:61–68`
It tests exact equality of column-major bytes. ProFit 2.7.2 normalised `Data$psf` in place at setup, so 60/66 accepted kernels differing from the cached `psf.rds` kernel is the expected signature of normalisation, and the 13 accepted-vs-source differences may be the same. **Add** a proportionality classification: `RENORMALISED` if `max|a/Σa − b/Σb| < 1e-10`, `RESIZED` if dimensions differ but the central overlap is proportional, else `DIFFERENT_KERNEL`; report counts. Until then the 13 are **UNRESOLVED**, not evidence of a defect.

### P2-2 Recovery-test tolerances too loose to catch regressions — `tests/test_synthetic_recovery.R`
Keep the test, add my original (`test_synthetic_recovery.R` from the V02 review: 0.08 mag, 8 % bulge Re, 5 % disc Re, 0.03 in q, 2°) as a second, stricter case at its own noise level.

### P2-3 Bundle manifest globs — `scripts/hectorbd_make_claude_review_bundle.py`
`outputs/*spatial*/*.tsv` is one directory level too shallow (the TSVs are `outputs/step2_spatial_calibration_*/<id>/cells_*.tsv`); `outputs/*pilot*/*.tsv` would only ever match the file that P1-2 prevented. Harmless.

---

## 4. Findings by requested item

### Item 3 — PSF geometry: see P1-1. Decision: **defect in code; fix before any observed PSF is built.**

### Item 4 — Historical P records
VERIFIED_FROM_BYTES: 68/68 `finesample = 1`, `psftype = 'empirical'`, ProFit 2.7.2, in both the R53 accepted layer and the original production layer. **Consequence:** the ProFit point-source bug at finesample > 1 (reproduced by me in 2.7.2 and 2.8.2) **did not affect any historical P fit**; the v1 flag `P_FINESAMPLE_GT1` is empty. The 13 accepted-vs-source and 60 accepted-vs-cache kernel differences are UNRESOLVED pending P2-1; the most probable explanation is in-place normalisation. The “failed final extraction helper” is P1-3: a provenance mismatch, not a data problem.

### Item 5 — Rfits comparison: see P1-2. **Orientation verified; comparison bug.**

### Item 6 — Sky, noise, calibration (VERIFIED_FROM_BYTES unless marked)
- Pedestals −0.0002…−0.0009 for 10/20/40 px masks (−0.02 to −0.06 σ): the images are background-subtracted to well below the noise; a fixed measured pedestal is the right policy.
- Blank-sky residuals are Gaussian to the precision measured: std 1.007–1.016, skew ≤0.03, 3σ tail 0.36–0.43 % (Gaussian: 0.27 %). The 3σ excess is consistent with faint undetected sources in “sky” pixels.
- ACF, all three fields and all three mask radii (my recount from the `acf_d*_raw.npz` files; the lag grid is integer and centred, and both axes agree):

  | Field | ρ(1) x / y | ρ(2) | plateau, lags 20–64 | Σρ, lags ≤5 | Σρ minus plateau |
  |---|---|---|---|---|---|
  | `11690105411641` | −0.075 / −0.071 | +0.053 | 0.003 | 2.19 | 1.86 |
  | `41691060304604` | −0.078 / −0.074 | +0.055 | 0.006 | 2.73 | 2.03 |
  | `11690092907692` | +0.014 / +0.024 | +0.140 | 0.068 | 12.8 | 4.64 |

  (10 px mask; 20 and 40 px give the same pattern.) **Decomposition (VERIFIED_FROM_BYTES):** the first field equals the other two plus a constant ≈0.09 at every short lag (−0.075+0.09≈0.015; 0.053+0.09≈0.14; 0.006+0.09≈0.095). So there are two components: (i) a weak, stationary short-range pixel covariance common to the fields, negative at lag 1 and positive at lag 2 — the regular alternating pattern expected from interpolation in resampling (INFERENCE; the producer's kernel would confirm it, §8 Q4) — giving N_eff ≈ N/2 for model comparison; and (ii) a field-specific long-range plateau that is large-scale background structure, not pixel noise: it barely falls with a 40 px mask and is not a plane, so it is a sky-systematic issue for that field, to be handled by a local background term or flagged, never folded into N_eff. **Correction to my earlier reading:** I first read the first galaxy alone as N_eff ≈ N/3–N/4; that conflated the plateau with the pixel covariance. **Requirement before any BIC:** use N_eff from the short-range component on fields without a plateau (≈N/2 here, to be confirmed on the 20 pilot fields), and flag fields whose plateau exceeds ~0.02. The plane-fit coverage failure at 40 px is correctly UNMEASURABLE, and a plane is in any case the wrong model for this structure.
- Gain: 1162 / 784 / 406 are `profoundGainEst` lower bounds (`step2_calibration.py:112`, field named accordingly). The spread anticorrelates with galaxy brightness (brightest, 41691060304604, lowest): galaxy structure inflates the variance and biases the estimator low. **Discriminating measurement (not a guess):** residual-variance versus model-brightness in annuli after a model fit of each galaxy and of its PSF stars; slope = 1/g. This is only possible after the first fit, so gain remains **UNRESOLVED** for now and cannot be set `VERIFIED_FROM_BYTES`.
- Units/zeropoint: BUNIT, MAGZERO, GAIN, EXPTIME absent from all headers. Segment magnitudes at an assumed 23.9 are 17.8, 16.7, 14.7 — plausible for Hector targets but not proof. **UNRESOLVED** until the producer metadata (§8) or a like-for-like comparison with WAVES segment photometry (`mag_rt`) is available. Do not set the two fields `VERIFIED_FROM_BYTES` from plausibility.

### Item 7 — Pilot: **PASS** (table in §1). The eligibility threshold (WCS ≤ 0.01″) is stricter than needed for science but right for a clean pilot. The first-galaxy candidate `11690092907692` is an anchor and `11690105411641` (calibration galaxy) is in the stratified ten; both have zero WCS offset and identical segmentation.

### Item 8 — The seven missing items
1–2. `outputs/*spatial*/*.tsv`, `outputs/*pilot*/*.tsv`: **harmless glob misses** (P2-3). 3. `step3…/summary.json`: **genuinely missing**, caused by P1-3. 4–7. `pilot_inputs_1453635/{originals_after.json, candidate_path_index.tsv, candidate_step0.json, SHA256SUMS.txt}`: **genuinely missing**, caused by P1-2. Nothing else is absent.

---

## 5. Gate decisions (original plan `CLAUDE_BD_RESTART_01` §14; E-items §7)

| Gate | Decision | Basis |
|---|---|---|
| G0 Step 0 | **BLOCKED** | science_image / original_segmentation / wcs_alignment / catalogue_label_mapping are satisfiable from existing files (Step 1 rows, catalogue job 1427444); units_zeropoint, gain_calibration UNRESOLVED (§4 item 6); upstream_sky and noise_protocol satisfiable from job 1437317 (pedestal, σ_MAD, z-statistics) |
| G1 Phase A lineage | **PASS** for the question it answers | 832/832 saved-vs-delivered segmentation; 822 identical, 10 pruned ≤4 px; likelihood support follows deterministically from `region = segim==label & !mask` (ProFuse `2a88cb7`), verified for the 3 dry-run galaxies by `Data$region` |
| G2 | **NOT_RUN** on observed data | code verified on synthetic |
| G3 | **PASS** (synthetic) | provenance tests at `98f6e08`+ |
| G4 F1–F9 | **PASS** | job 1443677 |
| G5 PSF E1–E5 | **FAIL → fix P1-1, then rebuild** | E-1, E-2 in code (`psf_v14.R`); E-3 encircled energy recorded by the builder (VERIFIED_FROM_CODE); E-4 per-object record exists in format; E-5 F6 PASS; but every rendered kernel is off-centre until P1-1 is applied. E-6/E-7 NOT_RUN |
| G6 models | **PASS** (declared) | `R/models.R`: S 7, P 7, D_FIXED4 8, BD_FREE_N 11 |
| G7, G8 | **NOT_RUN** | — |

---

## 6. Go / no-go

- **First observed fit of `11690092907692`: NO-GO today.** Required first: P1-1 applied and `test_psf_centering.R` passing on Friday; the galaxy's PSF v14 rebuilt and its row written to `config/psf_g4_pilot.tsv`; P1-2 applied and the pilot-inputs job rerun to produce the verified index row and `candidate_step0.json`; Step 0 evidence filled for the six fields that existing files support, with `units_zeropoint` and `gain_calibration` **remaining REPORTED** — which means the code gate will still refuse. That is correct behaviour. The honest route to a first fit is therefore either (a) producer metadata that certifies units and gain (§8), or (b) an explicit, documented *provisional* Step-0 policy for those two fields (e.g. `PROVISIONAL_ASSUMED_PENDING_PRODUCER`), implemented as a distinct status that `bd_step0_gate()` accepts only when `run_pilot.yaml` carries `provisional_calibration_acknowledged: true` and the galaxy count is ≤ 1. I recommend (b) only if (a) cannot be obtained within the week, and only with that code change plus a test showing the full-832 gate still refuses it.
- **20-object pilot: NO-GO** until the first fit has been inspected and the PSF rows exist for all 20.
- **832 run: NO-GO.**

---

## 7. Shortest next sequence (all SLURM, read-only except the two patches)

1. **Apply P1-1 and P1-2** (`patch -p1` from the repo root), commit, then run `slurm/psf_geometry.sbatch` and the pilot-inputs job again. Acceptance: all 21 geometry cases centred; `pilot_inputs_*/SHA256SUMS.txt` exists. Cost: minutes.
2. **Rerun `slurm/step3_historical_p_review.sbatch` from the committed tree** (P1-3) with P2-1's proportionality classification if ready. Acceptance: `summary.json` present; helper SHA in `binding_evidence.json` equals the commit; the 13 kernel differences classified. Cost: ~10 min.
3. **Fill Step 0 from existing evidence** for the six supportable fields (hashes of the Step 1 catalogue row, the WCS row, and `calibration.json`), leave two fields REPORTED, and write the producer question list (§8). Then decide route (a) or (b) above.

No new 832-scale job is needed; the F2 job must not be rerun.

---

## 8. Questions for the WAVES tile producers (only what certification needs)
1. What are the pixel units of the 0.30″ r-band cutouts (μJy, counts, counts s⁻¹, or scaled), and what AB zeropoint applies (is 23.9 correct)?
2. What sky treatment was applied before cutting (ProFound sky map? what box/grid), so that the residual pedestal we measure (≈ −0.0003) is interpreted correctly?
3. Is there an exposure/weight map for these coadds, and how should the per-pixel Poisson term be formed — i.e. the effective gain in (counts per image unit), or the number of coadded exposures and per-exposure gain?
4. What resampling kernel produced the 0.30″ grid from the 0.20″ native pixels (we measure ρ(2) > ρ(1) in the noise autocorrelation), so that the pixel-noise covariance can be modelled?
5. Which segmentation product is delivered (detection `segim_orig` or dilated total-photometry `segim`), so that segment-summed fluxes can be compared like-for-like with `mag_rt`/`mag_rc`?

---

## 9. Evidence ledger
- Repository: `489e73742b9cfe0e303c4c4162cf48bf5e74526f` (commit "Publish calibration and pilot preparation sources for scientific review", 2026-11-10T22:38+10:00); `SHA256SUMS.txt` PASS; `PUSHED_SOURCE/` identical.
- Jobs read: 1443186, 1443187, 1443677 (runtime/F1–F9); 1453593 (synthetic); 1453626 (PSF geometry); 1453616 (historical P); 1453629 (pilot); 1453634/35 (pilot inputs); 1437317 (calibration); 1453614/15/39/40 (spatial); 1427058/1427370/1427444 (Step 1, earlier).
- My commands: `Rscript /tmp/psf_centre_check.R <tree>` (centring, both trees); `Rscript tests/test_psf_v14_synthetic.R`, `tests/test_psf_branch_contract.R` on the patched tree; `python3 test_rfits_sample_tolerance.py`; `diff -r` clone vs bundle; JSON/TSV recounts above.
- Deliverables: `claude_psf_centre_fix.patch`, `claude_rfits_sample_tolerance.patch`, `test_rfits_sample_tolerance.py`, this report; SHA-256 in `SHA256SUMS.txt`.
