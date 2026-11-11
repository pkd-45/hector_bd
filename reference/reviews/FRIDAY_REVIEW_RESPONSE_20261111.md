# Response to the Claude review of 489e737 — Friday, 2026-11-11

Repository `pkd-45/hector_bd`, branch `main`. Starting commit
`489e73742b9cfe0e303c4c4162cf48bf5e74526f`. Review archive:
`reference/reviews/claude_20261111/` (all four entries of its `SHA256SUMS.txt`
verified; the uploaded zip held the same five files). This record does not
authorize observed fitting, the 20-object pilot or the 832 run.

## 1. What was applied

| Review item | Action | Source commit |
|---|---|---|
| P1-1 PSF centre | `claude_psf_centre_fix.patch` applied unchanged (clean, no fuzz). A comment at the two star-stamp sites records why `n/2` is right there (§3.2). | `2661a71` |
| P1-2 Rfits samples | `claude_rfits_sample_tolerance.patch` applied, then the comparison moved into `check_rfits_samples()` with the same tolerance (rel 1e-12, abs 1e-15). Added: sample-count check (`zip` truncated a short R list silently) and finiteness-flag check. | `2661a71` |
| P1-2 regression | Tests in `tests/test_prepare_pilot_inputs.py` call the production function `run()` uses: real 15-digit values accepted, transposed read rejected, coordinate/count/finiteness mismatches rejected, and the actual job-1453635 JSON evidence passes while the old `!=` rule fails on it. Against the unpatched `489e737` script, 8 of the new tests fail. The supplied `test_rfits_sample_tolerance.py` (which re-implements the rule) is kept unchanged in the archive and passes. | `2661a71` |
| P1-3 helper provenance | `step3_historical_p_review.py` refuses to run unless it and the R helper equal their committed blobs, and records `git_HEAD`, blob ids and tracked changes. | `2661a71` |
| P2-1 kernel classes | Helper adds value-level classification (`EXACT`, `RENORMALISED`, `RESIZED`, `DIFFERENT_KERNEL`, `DIFFERENT_KERNEL_AND_DIMENSIONS`, tolerance 1e-10 on normalised kernels) and physical centroid geometry for every accepted/source kernel. The original hash comparison is unchanged. | `2661a71` |
| Star-stamp convention | New `tests/test_psf_star_stamp_centre.R`, run by `slurm/psf_geometry.sbatch` with `test_psf_v14_synthetic.R` and `test_psf_branch_contract.R`. | `2661a71`, `817955c` |
| P2-2, P2-3 | Not done (P2-2 needs a file not on Friday, §6; P2-3 is a harmless bundle-glob issue). | — |

## 2. SLURM jobs (all from a clean committed tree)

| Job | Purpose | Commit | Result |
|---|---|---|---|
| 1458866 | PSF geometry + star-stamp | `2661a71` | `test_psf_centering.R` PASS; star-stamp test FAILED on combined assertion A4 (§3.3); later tests not reached |
| 1458867 | Pilot-input unit tests | `2661a71` | 12/12 PASS, none skipped |
| 1458868 | Pilot inputs, W11690092907692 | `2661a71` | COMPLETED; all outputs present, `SHA256SUMS.txt` verifies |
| 1458869 | Historical P metadata (68) | `2661a71` | COMPLETED in 2 min 11 s; `summary.json` present |
| 1458949 | PSF geometry + star-stamp + v14 helpers + P branch | `817955c` | All PASS: centring 21/21, star-stamp A1–A5, `D5 PSF synthetic helpers PASS`, `PSF_BRANCH_CONTRACT_PASS`; 62 min (submitted with `--time=02:00:00`; the script's 1 h default is too short) |

`R/psf_v14.R` is the same blob (`ead4ca5`, SHA-256 `1e4b22f2…`) at both commits.
F1–F9, F2 and the synthetic recovery suite were not rerun: none of them calls a
`psf_v14` render function.

## 3. Results

### 3.1 PSF geometry (P1-1)
21/21 cases centred and 180°-symmetric (max |Δ| 3.6×10⁻¹⁵ px, asymmetry
4.9×10⁻¹⁶, |Σ−1| = 0). Before the patch (job 1453626): 6/21, Δ up to 0.50 px,
asymmetry 0.70.

### 3.2 Star-stamp convention (patch hunks at `.fit_one_star_gama_v11/v12`)
The patch description says the per-star hunks fix "the 1 px acceptance window".
They do not touch it; they move the initial guess and ±3 px bounds. Whether
`n/2` is right depends on the convention of `xc` when the stamp is cut at
`round(xc)`:

- ProFound 1.35.0 (`1d56e38`) `xcen` is physical (`meanwt(x − 0.5)`); without v9 the star sits at `(n+1)/2 ± 0.5` in the stamp.
- `profuse_psf_v9.R` (SHA-256 `af6b8f28…`, sourced first by `build_pilot_psf_candidate.R`) defines `.refine_centroid`, which returns 1-based index coordinates; in production the star therefore sits at `n/2 ± 0.5`.

Synthetic Moffat stars (FWHM 3.2 px, con 2.7, total flux 10⁴, pixel noise σ = 1)
at seven sub-pixel positions, both paths, both fitters (job 1458949):

- A1 ProFound centre equals the physical truth within 0.05 px — PASS.
- A2 with v9, the refined centre equals truth + 0.5 (index convention) within 0.1 px — PASS.
- A3 without v9, the ProFound centre is used unchanged — PASS.
- A4 the converged v11 fit recovers the true centre within 0.01 px (max 0.0030) in both paths, although the initial guess was 0–1 px from the star — PASS. The patched initial guess does not bias the fitted centre.
- A5 production path: the star lies −0.45…+0.45 px from the patched `n/2` (it was −0.95…−0.05 px from the old `(n+1)/2`) — PASS. The patched guess and ±3 px bounds are centred on the expected star.

**Suitability window (not changed, needs a decision).** `.gama_suitability`
and `.gama_suitability_dr4` (`R/psf_v14.R:962`, `:1030`) still test
`|xcen − (stamp_size+1)/2| ≤ 1 px`. In production the fitted centres fall at
−0.95…−0.05 px from that reference: all 14 production fits pass, but a star whose
refined centroid lies near −0.5 px within its pixel is only 0.05 px from the edge,
while the other side has 1.05–1.95 px of margin. The window spans −0.5…+1.5 px
about the expected centre, so real-noise stars will be rejected preferentially on
one side. Centring it on `stamp_size/2` would match the patch, but it is a star
selection criterion and is left for the reviewer.

### 3.3 Star-fit optimiser (new finding, not caused by the patch)
`build_psf_for_image_v13` defaults to `profit_optim_dr4 = "L-BFGS-B"`, i.e.
`.fit_one_star_gama_v12` (`factr = 1e7`, `pgtol = 1e-7`, `R/psf_v14.R:888`).
Paired fits of the identical stamps above (14 per fitter, job 1458949):

| | v11 (BFGS) | v12 (L-BFGS-B, production) | truth |
|---|---|---|---|
| FWHM (px) | 3.193–3.208 | 3.240–3.309 (+1.3 to +3.4 %) | 3.2 |
| con | 2.668–2.711 | 2.88–3.55 (+0.18 to +0.85) | 2.7 |
| axrat | 0.998–1.000 | 0.968–0.982 | 1 |
| reduced χ² | 0.30–0.68 | 1.09–8.06 (higher than v11 in 14/14) | — |
| max centre error (px) | 0.003 | 0.032 | — |
| `.gama_suitability_dr4` pass | 14/14 | 7/14 (all 7 failures on χ² ≥ 3) | — |

All 28 fits report convergence code 0 and no parameter at a limit. The same inputs
on two nodes (job 1458866 on senv-c01, job 1458949 on smp-0-11) give v11
endpoints that agree to 2×10⁻⁵ px, but v12 endpoints that differ by up to 0.021 px
and 0.85 in reduced χ². A "converged" endpoint that moves with floating-point
details is not at the minimum.

History: job 1458866 asserted that both fitters recover the centre within 0.05 px;
v12 reached 0.051 px and the job failed. Commit `817955c` split that assertion:
A4 now tests the converged v11 fit, and v12 convergence is recorded in
`star_stamp_centre_report.json → v12_convergence_finding` rather than asserted
or hidden by a looser threshold.

Consequence (INFERENCE): v12 star fits are biased towards a broader, more
concentrated, slightly elliptical Moffat, which would propagate into the
median-parameter PSF, and well-measured stars can fail the χ² suitability test.
The cause is not isolated (stopping criteria, the bounded logit axrat near its
0.999 cap, or the bounded `bg`/`mag` are candidates), and the size of the effect at
realistic stellar S/N is UNRESOLVED; the high S/N here amplifies χ². The optimiser
choice and its settings are PSF-method decisions and were not changed. The test
gives a ready check for any change (e.g. `profit_optim_dr4 = "BFGS"`, or tighter
v12 stopping).

### 3.4 Pilot inputs (P1-2)
Job 1458868 wrote `originals_after.json` (`originals_unchanged: true`),
`candidate_path_index.tsv`, `candidate_step0.json` (status `BLOCKED`, all
calibration fields `UNRESOLVED`, `observed_fitting_authorized: false`) and
`SHA256SUMS.txt`. Copy SHA-256, numeric array hashes, originals and the Rfits
products are identical to the failed job 1453635; only the comparison changed.
The executed script hash equals the `2661a71` blob.

### 3.5 Historical P (P1-3, P2-1)
- Provenance: `git_HEAD 2661a71`, no tracked changes, driver and helper blobs equal HEAD.
- 68/68 accepted fits: `Data$finesample = 1`, `Data$psftype = "empirical"` (unchanged from 1453616).
- Accepted vs source kernel: 55 `EXACT`, **13 `RENORMALISED`**. In all 13 the dimensions are equal, both sums are 1 to 15 digits, max raw difference ≤ 4.2×10⁻¹⁷, max normalised difference ≤ 1.4×10⁻¹⁷ (≤ 1.7×10⁻¹⁶ of the peak), spread over 78–80 % of pixels: last-bit differences from renormalising a kernel whose sum is 1 ± ULP. They are the same kernels; no defect.
- Accepted vs cached `psf`/`psf_gama`: 60 `RENORMALISED`, 6 `EXACT`. Cached `psf_default` (the fallback, 33×33): 61 `DIFFERENT_KERNEL_AND_DIMENSIONS`, 5 `DIFFERENT_KERNEL` — expected, it was not the kernel used.
- **All 68 accepted historical kernels are off-centre by +0.485…+0.498 px in x and y.** This confirms from the bytes the review's inference that historical v13 PSFs carried the P1-1 offset (slightly below 0.5 because the Moffat wings are truncated).

## 4. Gates and data
No change to `config/step0.yaml`, `config/run_pilot.yaml`, `config/psf_g4_pilot.tsv`,
approval flags, suitability criteria, model definitions or tolerances. All writes
went to new `outputs/*_<job>` directories. Original compressed FITS files were
hashed before and after (unchanged); historical R53/R61 RDS files were read only.

## 5. Outstanding scientific blockers
1. Step 0: `units_zeropoint` and `gain_calibration` remain `REPORTED`/`UNRESOLVED`; `bd_step0_gate()` correctly refuses. Needs producer metadata (review §8) or the documented provisional route (review §6 b), which is a reviewer decision.
2. PSF: the target's v14 PSF has not been built and `config/psf_g4_pilot.tsv` is empty. Before building it, decide: (a) the suitability centre window (§3.2), and (b) the v12 convergence and shape bias (§3.3), since v12 is the production star fitter. `slurm/pilot_psf_candidate.sbatch` can then run on `outputs/pilot_inputs_1458868`.
3. E-6/E-7 and the original E1–E7/G5 criteria remain NOT_RUN.
4. N_eff / sky plateau policy (review §4 item 6) is not yet implemented.
5. P2-2 stricter recovery test is not added (file not on Friday).

## 6. Files needed from outside Friday
- Claude's original V02-review `test_synthetic_recovery.R` (0.08 mag, 8 % bulge Re, 5 % disc Re, 0.03 q, 2°) for P2-2. The README records it was never attached.
- Producer answers or metadata for review §8 (units/zeropoint, sky treatment, weight/exposure maps, resampling kernel, segmentation product).
- Optional: the reviewer's `/tmp/psf_centre_check.R` (not needed; `tests/test_psf_centering.R` covers it).
