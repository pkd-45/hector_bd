# Hector/WAVES bulge–disc restart — where everything is

**Prepared by:** Claude, 30 October 2026.
Status tags: **[V-today]** I opened it today · **[V-earlier]** I opened it in an earlier review in this conversation · **[R-ChatGPT]** listed in ChatGPT's file 02, not verified by me · **[R-provenance]** a path printed inside a bundle I reviewed, not opened · **[NOT REACHED]** I tried and could not access it.

I have no access to your Mac, to Friday, or to ChatGPT's library. Every path on those systems is historical evidence, not a confirmed present-day location.

---

## A. Upstream code — exact versions [V-today]

All repositories were queried live today. ProFit, ProFound, ProFuse and Highlander were cloned with full history.

| Repository | Current version | HEAD commit (date) | Historical runtime |
|---|---|---|---|
| github.com/ICRAR/ProFit (official) | 2.8.2 | `2dca6eb` (26 Oct 2026) | **2.7.2 = `a804a69`** (7 Jul 2026) |
| github.com/asgr/ProFit | 2.6.2 | — | Old fork; **do not use** |
| github.com/asgr/ProFound | 1.35.0 | `1d56e38` (26 Oct 2026) | Not recorded — needs a runtime fingerprint |
| github.com/asgr/ProFuse | 0.6.4 | `2a88cb7` (7 Jul 2026) | 0.6.4 = `2a88cb7` (same code) |
| github.com/asgr/Highlander | 0.4.0 | `a44df8d` (11 Sep 2026) | 0.3.4 (in the repository from 11 May to 16 Jun 2026) |
| asgr/ProPane | 1.12.3 | — | — |
| asgr/ProSpect | 1.13.3 | — | — |
| asgr/ProTools | 0.5.0 | — | — |
| asgr/ProUtils | 0.1.6 | — | — |
| asgr/Rfits | 1.24.1 | `9c951ff` (7 Oct 2026) | — |
| asgr/Rwcs | 1.9.1 | — | — |
| asgr/ParmOff | 1.0.3 | — | **Now imported by ProFit 2.8.2** |
| asgr/magicaxis, celestial, LaplacesDemon, cmaeshpc, ProGeny | 2.7.0, 1.5.7, 16.1.8, 1.0.1, 0.8.3 | — | — |
| asgr/ProBot | 0.9.1 | `0ab3419` (2 Oct 2026) | — |

ProFit commits after your historical build: `e014094` (7 Oct, FFT speed and OpenMP), `04d9750` (version bump), **`a3b3262` (23 Oct, PSF in-place normalisation fix)**, `2dca6eb` (26 Oct, speed and memory work plus non-square PSF/pixel-scale fixes; no effect for square inputs).

Files that define the behaviour you care about (line numbers at the commits named):

| File | Commit | Lines | Content |
|---|---|---|---|
| ProFuse `R/profuseSingleImage.R` | `2a88cb7` | 1–41 | `profuseFound2Fit` defaults: `tightcrop=TRUE`, `autoclip=TRUE`, `roughpedestal=TRUE` |
| ″ | ″ | 42–46 | autoclip → NA |
| ″ | ″ | 62–98 | optional `cutbox` cutout |
| ″ | ″ | 124–150 | segim `NULL` → re-segment; supplied → `static_photom` |
| ″ | ″ | 176–202 | target label at `loc`; region = target label(s) minus mask |
| ″ | ″ | 205–215 | `tightcrop` bounding box |
| ProFound `R/profoundProFound.R` | `1d56e38` | 152–158 | static mode: `sky=0`, `iters=0`, `redosky=FALSE`, `redosegim=FALSE` (unchanged since `b5acc15`, 2021) |
| ProFit `R/profitSetupData.R` | `a804a69` | 230–231 | default region = label at image centre |
| ProFit `src/libprofit/src/profit/model.h` | `a804a69` → `a3b3262` | `set_psf` | in-place PSF normalisation → copy |

Vignette sources that mirror the RPubs posts (read today, narrative prose):

| Path | Relevance |
|---|---|
| `ProFit/vignettes/ProFit-PSF-Convolution-Can-Be-Convoluted.Rmd` (Taranu) | Discretised-kernel bias; finesample² scaling |
| `ProFit/inst/vignettes-nonCRAN/ProFit-Sigma-Maps.Rmd` | σ construction; two meanings of "gain" |
| `ProFound/vignettes/ProFound-Segim-Maps.Rmd` | Expanding and merging segments — operations we must avoid |
| `ProFound/vignettes/ProFound-Sky_Shenanigans.Rmd` | Sky algorithm; skyRMS = 50th − 15.9th percentile |
| `ProFuse/vignettes/ProFuse-Found-2-Fit.Rmd`, `ProFuse-MF2F-Setups.Rmd` | Found2Fit usage |

Not yet read: ProFound Blender, Colour, Complex-Fit, MultiBand-Stack, TipsTricks, Source-Finding, Beginners; ProFuse AGN-PSF and Multi-Images; ProFit Galaxy-Fit, Make-Model, PSF-Fit, Full-Monty, Isophotal; Highlander Quick-Fit and LD-Algorithms.

## B. Literature — what I actually read

| Paper | Identifier | My access today |
|---|---|---|
| Casura et al. 2022, GAMA bulge–disc on KiDS | MNRAS 516, 942; arXiv:2208.07608 | **[V-today]** §§1–4.3 in full (data, ProFound, segmentation, sky, σ, PSF, fitting, flags, model selection, `_SEGRAD` truncation) |
| Barsanti et al. 2021, bulge and disc colours in clusters | ApJ 911, 21; arXiv:2012.12480 | **[V-today]** §§I–V in full (cutouts, neighbours, iterative re-segmentation, PSF, Bayes-factor selection, filters, internal tests) |
| Barsanti et al. 2021, SAMI bulge/disc stellar populations (co-author S. M. Sweet) | ApJ 906, 100; arXiv:2011.04873 | Abstract and metadata only |
| Barsanti et al. 2022, spin–filament flips | arXiv:2208.10767 | **Not read** |
| Robotham et al. 2017, ProFit | MNRAS 466, 1513; arXiv:1611.08586 | Excerpts (point-source PSF handling; repository locations) |
| Robotham et al. 2018, ProFound | MNRAS 476, 3137; arXiv:1802.00937 | Abstract only; behaviour checked in source code instead |
| Robotham, Bellstedt & Driver 2022, ProFuse | MNRAS; arXiv:2204.04897 | Abstract; behaviour checked in source code instead |
| Bellstedt et al. 2024, ProFuse star-formation histories | MNRAS 528, 5452; arXiv:2307.02788 | **Not read** (the historical `profuse_select.R` cites its §3.3.1) |
| Casura PhD thesis 2022 | arXiv:2303.10077; ediss.sub.uni-hamburg.de/handle/ediss/9818 | **Not read** |
| Robotham et al. 2024, ProPane | arXiv:2311.01761 | **Not read** |
| Häußler et al. 2022, Galapagos-2/GalfitM/GAMA | arXiv:2204.05907 | **Not read** |
| Kaur et al. 2025, WAVES-Wide target selection | arXiv:2510.11132 | Excerpt: rizY detection; "total" (dilated) versus "colour" (detection) segments |
| Cook et al. 2024, WAVES star–galaxy separation | MNRAS 535, 2129 | Excerpt: ProFound segments; large fragmented sources manually regrouped |
| Bellstedt et al. 2020, GAMA KiDS/VIKING ProFound photometry | arXiv:2005.11215 | Excerpt only; the WAVES pipeline is described as "similar" |
| Bellstedt et al. (in prep.), WAVES photometry | — | **Unpublished**; needed for S0.1–S0.4 — ask the WAVES team |
| Cook et al. 2025, DEVILS size–mass of components | MNRAS (2025), doi:10.1093/mnras/staf575 | **[NOT REACHED]** abstract only |

## C. RPubs (Aaron Robotham) — https://rpubs.com/asgr

Index **[V-today]**: 41 posts. RPubs retires in June 2027; documents stay readable until 31 Dec 2031, so save local copies. **Post bodies are [NOT REACHED]** — they load in an embedded frame my fetch tool cannot open. Where a matching vignette exists in the repositories, I read that source instead (§A).

| ID | Title | Relevance |
|---|---|---|
| 274425 | ProFit – Sigma Maps | High (σ; source read) |
| 296847 | ProFound – Sky Shenanigans | High (source read) |
| 296848 | ProFound – Source Finding | Medium |
| 333281 | ProFit – PSF Convolution Can Be Convoluted | **High** (source read) |
| 385865 | Prior Knowledge | High for Bayes-factor selection |
| 715082 | ProFit: Found 2 Fit | High (source read) |
| 844057 | ProFuse: MF2F Setup | High (source read) |
| 806303 | ProFuse: Simple Fit | Medium |
| 1015573 | ProFound: Segmentation Maps | High (source read) |
| 1015572 | ProFound/ProFit: A Complex Fit | Medium |
| 1015936 | ProFit: Galaxy Fitting Example | Medium |
| 1015938 | ProFit: PSF Fitting Example | Medium |
| 1015578 | Lost and ProFound: Tips and Tricks | Medium |
| 1015569 | ProFound: Colour Me Happy | Medium (`segim` versus `segim_orig`) |
| 1367546 | ProFound: It's All in the Blend | Medium (deblending) |
| 307935 | ProFit/ProFound – The Full Monty | Medium |
| 286460 | ProFit/ProFound – Using Isophotal Ellipses | Low–medium |
| 727397, 718358 | ProFit: Multi Image Mode; Offset Images | Later extension |
| 1427940 | MCMC Algorithms in LaplacesDemon | Inference |
| 1433553 | ParmOff | ProFit 2.8.2 dependency |
| 1328227, 1328228, 1422211, 1012748 | Stacking and warping (ProFound, ProPane) | Noise covariance |
| 1376433 | ProUtils: Don't be Noisy | Stacking noise |
| 1399576, 1139429, 274701, 274942, 274426, 275014 | Beginner's guide, Rfits, Make Model, README, OpenCL, SAMI cluster demo | Low |
| Others | ProSpect, ProGeny, ProBot, Arrow, EBL | Not relevant |

## D. Project data and code on Friday — historical paths

| Path | Status | Role |
|---|---|---|
| `/data/scratch/uqpdas/WAVES_bd_v32_production842/baseline/FIXED4_seed666/<uberID>/r_1333_v25/` | [R-provenance] (printed in v04 card provenance; ChatGPT gives one exact example) | Original production: `inputs.rds`, `fit_FS.rds`, `fit_PD.rds`, `fit_BD.rds`, `psf.rds` |
| `…/baseline/GAMA_seed666/` | [R-ChatGPT] | Free-n branch |
| `/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R53R1_FULL832_PRODUCTION_EXECUTION_PATH_REPAIR_20260921/runs/<uberID>/<model>/<start>/accepted_fit.rds` | [R-provenance] root; [R-ChatGPT] template | Accepted R53 fits behind the catalogue |
| `/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R54R1_…_20260922/results/R54_BEST_MODELLIST_LONG.tsv` | [R-ChatGPT] | Best-endpoint model list |
| `/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_…_20260930/results/R61_NATIVE_MODEL_SELECTION_832.tsv` | [R-ChatGPT] | Native S/P/D_FIXED4 selection |
| `/data/uqpdas/WAVES_v32/4_run_production_v32.R` | [R-ChatGPT]; **source not recovered** | Production runner — **highest-priority file** |
| `/data/uqpdas/WAVES_v32/profuse_bdmodels_v28.R`, `…_v19.R`, `select_v32_primary_fixed4.py` | [R-ChatGPT]; not recovered | Model wrappers and selector |
| `/data/uqpdas/WAVES_v38_profit272_pilot/Rlib` | [R-ChatGPT] | Pinned historical R library — fingerprint the versions here |
| `/home/uqpdas/.conda/envs/waves_r/bin/Rscript` | [R-ChatGPT] | Historical R |
| `/data/scratch/uqpdas/g011q1_production842_static/canonical_inputs` | [R-ChatGPT] | Canonical image/σ/mask/PSF used by R67/R74 (matched-engine branch) |
| `/scratch/uqpdas/WAVES_BD_PROFIT_GALFIT_C0R4C16R74R3R1_REAL832_DUAL_ENGINE_BDFN_EXECUTION_20260929` | [R-ChatGPT] | Dual-engine BD_FREE_N run |
| `/data/uqpdas/WAVES_BD_CROSSCODE_V010_TOOLS/bin/galfitm-1.4.4-linux-x86_64` | [R-ChatGPT]; SHA `bdfc8d20…` [V-earlier] | GALFITM binary |

## E. Mac paths — historical

| Path | Status | Role |
|---|---|---|
| `/Volumes/waves/WAVES/<uberID>/<uberID>_r_1333.fits.gz` | [R-ChatGPT] | **Source images (1333×1333)** — and possibly the segmap. Confirm via `inputs.rds$segmap_path` |
| `/Volumes/waves/WAVES_bd_v32_production842/` | [R-ChatGPT] | Mirror of the production tree |
| `~/Downloads/Hector_input_catalogue_v03.2.csv` | [R-ChatGPT] | Hector input catalogue |
| `~/Downloads/HectorBD_Rband_v01_FINAL_QC1_20261008T202048.tar.gz` | [R-ChatGPT] | Frozen QC1 release |
| `~/Downloads/HectorBD_MGE_OBSERVED_ELLIPSE_MATCHED_AUDIT_*` | [V-earlier] (named in the v03 rebuild script) | MGE ellipse audit |
| Conda environment `hector-photometry` | [V-earlier] (named in the v03 script) | Python environment for the audits |

## F. ChatGPT library

ChatGPT's file 02 lists 221 library items with SHA-256 values. **[R-ChatGPT]** — I cannot open that library. The most relevant items:

- the frozen QC1 tarball (`b3633cc8…`)
- the September raw-evidence zip, which contains `profuse_psf_v13.R` and `profuse_select.R`
- `WAVES_BD_team_explainer.tex`, which describes the island-retention step
- the July master handoff DOCX, which contains the smoke72 helper text

## G. Claude review outputs in this conversation [V-today]

These are under `/mnt/user-data/outputs/` (download them; this location is not permanent):

- the C2′ to C2⁷ audits (8–9 Sep)
- the science reviews: S_SINGLE/BD_FIXED4 (18 Sep), R23–R38 (21 Sep), R64 (24 Sep), R68 (25 Sep), R69 (30 Sep), R69R3R1 (1 Oct)
- the three-way comparison review and its figures (`HectorBD_review_figures/`)
- the 13-survivor re-review and the V03 ellipse final check (6 Oct)
- `HectorBD_RE_RADIAL_review/CORRECTED_TRUNCATION_TABLE_39.csv`
- `r68_stage0_reanalysis.py`
- the presentation plans v1 and v2
- this restart set in `bulge-disk-restart/`

## H. What to collect from Friday for Phase A

Run read-only commands on Friday, then bring back the outputs:

```bash
# 1. Where are the inputs and what do they point to?
ROOT=/data/scratch/uqpdas/WAVES_bd_v32_production842/baseline/FIXED4_seed666
ls "$ROOT" | head; ls "$ROOT"/*/r_1333_v25/inputs.rds | wc -l
/home/uqpdas/.conda/envs/waves_r/bin/Rscript -e '
x <- readRDS(Sys.glob("'"$ROOT"'/*/r_1333_v25/inputs.rds")[1]); str(x, max.level=1)
cat("image_path:", x$image_path, "\nsegmap_path:", x$segmap_path, "\n")'

# 2. What is in the segmap file? (HDUs, header, which product)
python -c "from astropy.io import fits; import sys; f=fits.open(sys.argv[1]); f.info(); print(repr(f[0].header))" <segmap_path>

# 3. Runtime fingerprints
Rscript -e 'for(p in c("ProFit","ProFuse","ProFound","Highlander","LaplacesDemon")) cat(p, as.character(packageVersion(p, lib.loc="/data/uqpdas/WAVES_v38_profit272_pilot/Rlib")), "\n")'

# 4. The production runner and loader
ls -l /data/uqpdas/WAVES_v32/; grep -n "profuseFound2Fit\|segim\|cutbox\|tightcrop\|autoclip\|side_fit\|island\|label" /data/uqpdas/WAVES_v32/4_run_production_v32.R
```

The fourth command alone may answer whether production passed the WAVES segmap, used a `cutbox`, or left `tightcrop` and `autoclip` at their defaults.
