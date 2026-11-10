# Hector/WAVES bulge–disc restart — what to do and what to check

**Prepared by:** Claude, 30 October 2026. Use with files 00 (review of ChatGPT's packet), 02 (file locations) and 03 (code inventory).
**Status:** a plan only. No fitting, segmentation, PSF or sky estimation was run for this document.

Labels: **[V]** verified by me from primary source or code · **[R]** reported elsewhere, not verified by me · **[I]** inference to be tested.

---

## 0. The contract

The galaxy is fitted on **every usable pixel of its WAVES target segment, and only those**. The segmentation is never altered.

Definitions, used throughout:

- `L` is the WAVES label map, exactly as delivered. Its file, HDU and SHA-256 are recorded.
- `A_g` is the documented set of labels that belong to galaxy *g*. This may be more than one label where WAVES regrouped a fragmented source.
- `T_g = (L ∈ A_g)` is the target membership in the **original** pixel grid, including disconnected islands and holes.
- `Q` is usable data: finite image, finite and positive σ, and not flagged by a documented detector mask.
- **`F_g = T_g ∧ Q` is the likelihood support.** Every pixel in `T_g` but not in `F_g` is written to an exclusion ledger with its coordinate and reason.

**Prohibited operations on `L`, `T_g` or `F_g`:** re-detection, dilation, erosion, hole filling, relabelling, merging, island removal, trimming by fitted radius, signal-to-noise or brightness clipping, residual-based masking, autoclip, and any per-model or per-engine segmentation. If the catalogue position falls on sky, stop that galaxy; never fall back to the nearest label.

**The first and blocking gate:** `F_g` must equal the original `T_g ∧ Q`, compared **in original coordinates**. Equality between two engines, or between two files we wrote ourselves, does not count. This is the gate missing from the earlier matched-engine reviews.

---

## 1. What is established so far

| # | Fact | Evidence |
|---|---|---|
| E1 | Production used ProFuse 0.6.4 `profuseFound2Fit`, which is identical to today's repository code | [V] ProFuse `2a88cb7` is both HEAD and the recorded historical version |
| E2 | Its fit region = every pixel carrying the target label in the segim it is given, minus the mask (all islands kept) | [V] `R/profuseSingleImage.R` line 199 |
| E3 | `tightcrop = TRUE` (the default) crops image, σ, segim and mask to the region's exact bounding box | [V] lines 205–215. This explains the historical stamp sizes: W176980350403709's frame was 25×35 = its 875 "broad" pixels |
| E4 | `autoclip = TRUE` (the default) sets extreme pixels to NA **before** the region is built; any target pixels it hits silently leave the fit | [V] lines 42–46 (thresholds 100× the 0.1%/99.9% quantile offsets) |
| E5 | With `segim = NULL`, ProFound **re-segments the image** (`sky = 0`). With a supplied segim, static mode: no dilation, no re-segmentation, `sky = 0` | [V] `profuseSingleImage.R` lines 124–150; `profoundProFound.R` lines 152–158 |
| E6 | Archived `fit$profound$sky = 0` does **not** reveal which branch ran — both give zero | [V] from E5 |
| E7 | Low-level `profitSetupData` defaults the region to the label under the image centre | [V] lines 230–231 @`a804a69` |
| E8 | ProFit 2.7.2 normalised the caller's PSF matrix in place (drift; saved PSFs are post-normalisation). Fixed 23 Oct 2026 | [V] `a3b3262` diff |
| E9 | PSF code: double MAD scaling in two places (star clipping; per-star σ fallback); linear PA clipping; size cap with zeroed corners and a 21-px floor loses 1–5% of the flux at β ≈ 2.1 across KiDS seeing | [V] code, plus my calculation (file 00, G2–G4, corrected 31 Oct) |
| E10 | Original loader, island filter, sky/gain code and production runner are **not recovered** | [R] ChatGPT file 03 C12; agreed |
| E11 | The project's earlier finite-size validation showed historical Sérsic Re are inflated by extrapolation beyond tight segments — exactly as Casura+2022 §3.3.3 documents | [V] earlier reviews; Casura §3.3.3 |

**What is not yet known:** whether production passed the WAVES map or `NULL`; whether a `cutbox`/`side_fit` cut the segment before Found2Fit; whether the island filter removed pixels; whether autoclip removed any. **Phase A answers all four from existing files, without fitting.**

---

## 2. Step 0 — decisions before any code (blocking)

| ID | Decision | How to settle it |
|---|---|---|
| S0.1 | **Which WAVES product is `L`:** detection segmentation (`segim_orig`-like), dilated total-photometry segmentation, or the regrouped/group map | Inspect the FITS headers and HDUs at `inputs.rds$segmap_path`; ask Sabine Bellstedt or the WAVES team which product it is and how regrouping is encoded |
| S0.2 | **Target mapping `A_g`:** label(s) per galaxy, from catalogue ID → label, never from a position lookup alone | WAVES catalogue segID/groupID columns; record the evidence |
| S0.3 | **Science image:** which band and file, pixel scale, units, zeropoint, and whether it is already sky-subtracted. WAVES imaging is ProPane-warped, so noise is correlated | FITS header plus WAVES documentation (Bellstedt et al. in prep.) |
| S0.4 | **σ source:** weight or inverse-variance map, or reconstructed? Which gain applies (vignette: two different "gains" exist)? | Header, WAVES products; never inherit `gain = 687` unverified |
| S0.5 | **Sample:** frozen 832, keeping the 842 and the 10 exclusions as a separate record | Existing lists |

Nothing downstream starts until S0.1–S0.3 are written down with evidence.

---

## 3. Phase A — read-only lineage audit of the existing fits

Purpose: establish *exactly* what pixels were fitted historically, and why. **No fitting.**

For every galaxy, from `inputs.rds`, `fit_*.rds$Data`, the accepted R53 RDS and the original FITS:

| Check | Pass condition | What a failure means |
|---|---|---|
| A1. Original `T_g` rebuilt from `L` at `segmap_path` | Label(s) found; hash recorded | Mapping problem; stop |
| A2. Map `Data$region` back to original coordinates (requires the cutout offset) | Offset recovered exactly | Coordinate lineage lost |
| A3. `N_missing = |T_g \ region|`, `N_extra = |region \ T_g|` | Both 0, or every missing pixel explained in the ledger | Pixels were dropped or added |
| A4. Classify each missing pixel | Into: cutbox/`side_fit` crop · island filter · autoclip NA · mask · invalid σ | Identifies the historical "cut" |
| A5. Does `region` match **any** WAVES label set? | Yes | No ⇒ ProFound re-segmented the stamp (`segim = NULL` branch) |
| A6. Does the segment touch the edge of the 1333×1333 source image? | No | The source image itself truncates the segment; flag it |
| A7. `Data$psf`: kind (analytic Moffat / empirical / default), dimension `d_use`, sum, hash | Recorded per galaxy | Unknown PSF provenance |
| A8. σ provenance: skyRMS and gain values used, units | Recorded | — |
| A9. Package fingerprints from the runtime: ProFit, ProFuse, ProFound, Highlander, LaplacesDemon (version and commit) | Recorded | — |

**Deliverable:** one table, 832 rows, columns for A1–A9, plus overlay PNGs for the pilot list (§12). This single product settles whether, where and how the pipeline departed from the segments.

---

## 4. Phase B — build fresh inputs

Per galaxy:

1. **Likelihood image** = a cutout of the source image covering the full bounding box of `T_g`, **plus padding** of at least max(PSF half-size, 3×FWHM). The padding never enters the likelihood. It is not needed for convolution correctness — libprofit pads the drawing domain by PSF/2 internally (`model.cpp`, `analyze_expansion_requirements`) — but it makes the render domain explicit, keeps segment edges visible on cards, and gives the sky diagnostics room.
2. **Region** = `F_g` exactly, passed **explicitly**. Never rely on a centre default (E7).
3. **No `autoclip`, no `cutbox` smaller than the segment, and `tightcrop = FALSE` by design.** `tightcrop` itself drops no pixels and libprofit pads for convolution, so it is a transparency choice, not a correctness fix. If Found2Fit is used at all, call it with `segim = L` (the original map), `tightcrop = FALSE`, `autoclip = FALSE`, `fit_extra = FALSE` unless neighbours are deliberately fitted, and `loc` set to a pixel inside `T_g`. Then assert that `Data$region == F_g`. The cleaner route is direct low-level `profitSetupData` with explicit inputs.
4. **Neighbours:** their pixels lie outside `F_g` anyway. Their **light leaking into** `T_g` is a declared variant: simultaneous fit when overlap and brightness criteria are met (Barsanti/Meert: R100 overlap and isophotal area > 5% larger). The primary run is target-only.
5. **Sky-diagnostic product:** a separate, larger cutout with its own object mask. This mask is a separate file; it **never** modifies `L`.
6. Hash image, σ, region, PSF and `L` at load, before the objective, after the fit and at export. All five must be unchanged.

---

## 5. Phase C — sky

1. Establish whether the WAVES image is already background-subtracted, and by what method (S0.3).
2. Measure any residual pedestal on background pixels outside **all** segments, beyond an object mask that is a separate dilated product. Use block bootstrap for the uncertainty, because the pixels are correlated.
3. Candidate method, Casura+2022 §3.1.3: box-wise medians on an aggressively dilated object mask, interpolated with a bicubic spline. ProFound's `skyRMS` is the 50th minus 15.9th percentile — a lower-half robust width.
4. **The primary fit holds the sky fixed.** Declared variants: sky ± its measured 1σ, and a free constant sky with k counted. A small target-only segment may not constrain a free sky; report that.
5. Watch for over-subtraction around large galaxies, which steepens outer profiles (Owers+2019, quoted in Barsanti+2021 §III.4).

## 6. Phase D — σ (noise)

1. Units and gain from S0.4. The object Poisson term is σ² = σ_bg² + max(model, 0)/g_eff, with `g_eff` valid **in the image's units**.
2. Check: blank-sky normalised residuals ≈ N(0,1), tails, spatial variation.
3. **Correlated noise from ProPane warping:** measure the autocorrelation and quote an effective N for any BIC or AIC.
4. **Mandatory probe:** with the model and `F_g` held fixed, doubling σ must divide χ² by exactly 4.

## 7. Phase E — PSF

| Item | Action |
|---|---|
| E-1 Double MAD (clip, line 630; σ fallback, lines 864 and 1133 — all in `BULGE_DISK_03_CODE_REVIEW.md`) | Fix both to a single 1.4826 factor; re-run star selection; report which stars change |
| E-2 PA clip | Clip on a circular statistic (doubled angle) |
| E-3 Size cap (21-px floor) and zeroed corners | Read `d_use`, FWHM and β from each `psf.rds`; measure encircled energy of the adopted kernel versus the analytic Moffat (bound: 1–5% lost at β ≈ 2.1). Cap so that ≥ 99% of the flux is retained, or document the loss; never renormalise silently |
| E-4 Kind | Record analytic, empirical or default per galaxy; the same kernel goes to fitting, display and scoring |
| E-5 Library mutation (E8) | Use ProFit ≥ `a3b3262`, or pre-normalise in double precision; hash the PSF before and after repeated evaluations — they must be identical |
| E-6 Validation | Leave-one-star-out residuals; spread across reasonable PSFs; propagate that spread to n, Re and B/T for the pilot |
| E-7 Spatial | Stars nearest by quadrant (Casura recipe); check FWHM variation across the 400″ field |

## 8. Phase F — renderer and likelihood tests (before any fit)

| Test | Pass condition |
|---|---|
| F1 Padding convergence | The model on `F_g` changes by ≪ σ as external padding grows, and matches libprofit's internally padded render of the tight stamp |
| F2 **Finesample convergence** (Taranu vignette) | The objective and parameters are stable from finesample 1 → 3 → 5 for compact bulges and P models |
| F3 Angle convention | A rotated, flattened test source: ProFit `ang` = from the 2nd matrix axis, counter-clockwise, direction (−sin θ, cos θ). Any custom aperture or overlay code must reproduce ProFit's own render |
| F4 Centre convention | An off-centre source in a non-square image: pixel *i* centre = *i* − 0.5 in native coordinates; no transposes |
| F5 Re semantics | `re` is the semi-major half-light radius; exponential Rs = Re/1.67835 converted in **one** place only |
| F6 PSF immutability | Hash identical before and after 100 identical evaluations |
| F7 σ probe | 2σ ⇒ χ²/4 |
| F8 Objective recomputation | The saved model plus saved arrays reproduce the reported likelihood exactly |
| F9 Region honoured | Changing one pixel outside `F_g` leaves the likelihood unchanged; changing one inside changes it |

## 9. Phase G — models

Declare before fitting:

| Model | Components | Free parameters (k) |
|---|---|---|
| S | Sérsic | x, y, mag, Re, n, q, PA (7) |
| P | point source + exponential | x, y (tied), mag_ps, mag_d, Re_d, q_d, PA_d (7) |
| D_FIXED4 | n = 4 bulge (round or free q, **declare which**) + exponential | historical: 8 |
| BD_FREE_N | Sérsic bulge + exponential, free n_b, q_b, PA_b | 11 |

- Bounds: write them down. Reference ranges are Casura Table 2 and Barsanti Table 1. Re upper bounds must not allow solutions far beyond `T_g` without a flag.
- Handle bulge–disc swaps explicitly (Casura §3.2.6). Never let the "bulge" be larger than the disc without flagging it.

## 10. Phase H — inference

- Pin one Highlander and LaplacesDemon version (Highlander 0.3.4 versus 0.4.0 changed prior handling). State whether `Data$prior` is used.
- Use several starting structures, not just seeds. Keep the best admissible objective under a declared rule.
- Uncertainties come from posterior samples (Casura found MCMC errors underestimate true errors by 2–3×). Start-to-start spread is not an error bar.

## 11. Phase I — model selection

- Options: Casura's DIC cuts calibrated by visual inspection; Barsanti's ln BF > 60 plus physical filters; the historical common-renderer BIC (R61). Thresholds are **not** transferable between methods.
- Marginal likelihoods depend on prior widths ("Prior Knowledge" vignette, 385865). Report prior sensitivity for any Bayes-factor rule.
- Calibrate on a pre-declared, visually inspected subset; validate on a separate one. Keep "ambiguous" and "unfittable" as outcomes.

## 12. Phase J — outputs and flags

**Primary sizes and magnitudes follow Casura §3.3.3:** the intrinsic profile integrated to the **segment radius** (`*_SEGRAD`), alongside the infinite Sérsic values, with both definitions stated.

Flags:

- Casura's: model flux inside the segment < 20% (flag) or < 50% (caution); "irregular segment" (segment-aperture magnitude minus segment-radius magnitude > 0.3); parameter at a bound; position offset > 2″; χ²ν > 80; B/T < 0.001 or > 0.999.
- Ours: segment touches the source-image edge (A6); PSF kind and encircled-energy loss; exclusion-ledger count; centre-to-segment-edge distance divided by R50.

Cards show the **actual segment contour**, excluded pixels, data, model and residual, and the correct PSF kind.

## 13. Phase K — validation and pilot

**Pilot of 20, declared before any fit.** Ten known cases:

| Galaxy | Why it's in |
|---|---|
| W42710765002459 | Sat at the historical stamp edge |
| W221990044410325 | Galaxy larger than its stamp |
| W351690654304280 | Bright neighbour |
| W176980350403709 | Tiny segment; Re twice the frame |
| W12690079611187 | ChatGPT's coverage case |
| W41691060304604 | Clean round elliptical (control) |
| W350680203312123 | Edge-on disc |
| W18680557409852 | Barred galaxy |
| W49710254502112 | Diffuse contamination |
| W213991001205426 | Round, low-n, low surface brightness |

Plus ten drawn at random across strata of size, n and crowding.

Regression tests after the pilot:

- The observed-plane semi-major comparison with MGE and ProFound, as finite-light radii, should reproduce or improve on the earlier ≈ 0.985 for well-covered galaxies.
- Mocks with **frozen segments**, generated by an independent renderer.

## 14. Acceptance gates (for the Claude reviewer)

| Gate | Requirement | Blocking? |
|---|---|---|
| G0 | S0.1–S0.3 written down with evidence | Yes |
| G1 | Phase A table complete for all 832; every missing pixel classified | Yes |
| G2 | New inputs: `F_g == T_g ∧ Q` in original coordinates; ledger complete | Yes |
| G3 | Hashes of `L`, image, σ, region and PSF unchanged at load, objective, after fit and at export | Yes |
| G4 | F1–F9 all pass | Yes |
| G5 | PSF items E-1 to E-5 resolved; E-6 reported on the pilot | Yes for the pilot |
| G6 | Models, bounds, k, and priors-or-not declared before the pilot | Yes |
| G7 | Pilot cards inspected; no fit uses pixels outside `F_g` | Yes before the full run |
| G8 | Selection calibrated and validated on separate subsets | Before any science claim |

## 15. Red flags for the new chat

Stop and escalate if any of these appear:

- `autoclip = TRUE`, any `cutbox`, `segim = NULL`, or `fit_extra = TRUE` without a declared neighbour policy in a fitting call. (`tightcrop = TRUE` is a design deviation to log, not a correctness failure.)
- A region built from anything other than the original `L`.
- "Support equality" checked only between engines or between our own files.
- Any dilation, merge or relabel function applied to `L` (`profoundMakeSegimDilate`, `profoundMakeSegimExpand`, `profoundDilate`, tolerance changes).
- A display-expanded image (V05/V05R1) used as a science input.
- A PSF described as "empirical" without `output_kind` evidence.
- Mixing the R54, R61 and R74 selection layers.
- A full 832 run launched before G0–G7 pass.

## 16. Order of work

**Step 0 → Phase A (read-only) → decide what changed → Phases B–F (inputs and tests, zero-fit) → pilot (Phases G–J on 20) → review → full run → Phases I–K.**

Each arrow is a stopping point for review.
