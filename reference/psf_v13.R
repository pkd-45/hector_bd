## ============================================================================
##  files/profuse_psf_v13.R
##
##  GAMA DR4 BDInputsRv05 PSF logic (Casura+2022; table id=698).
##  Fork of profuse_psf_v12: `build_psf_for_image_v13` uses **strict**
##  `select_clean_stars_v10` only (star_frac defaults 0.04, gri / skycut=2 prose),
##  `.gama_suitability_dr4` (no bg in 2-sigma clip), DR4 PSFQUAL bits **without**
##  v12's adaptive-pool misuse of bit 128.  Optional DR4 PSFQUAL bit 128 via
##  `psfqual_skycut_segstats_bounds = c(lo, hi)` row-count heuristic.
##  `profit_optim_dr4`: GAMA uses BFGS; default `L-BFGS-B` for WAVES stability.
##
## ============================================================================
##
##  Original v11 description (still applies to the unchanged parts) ----------
##  v11 PSF builder — Casura+2022 / GAMA DR4 BDInputsRv05 PSF logic (ProFound+
##  ProFit Moffat+BFGS; arXiv:2208.07608, MNRAS 516, 942; schema
##  https://www.gama-survey.org/dr4/schema/table.php?id=698 )
##  PLUS stability upgrades from this WAVES pipeline:
##    - sub-pixel centroid refinement (.refine_centroid from v9) on the
##      science image before the PSF stamp is extracted;
##    - segment-dilated local sky + ring RMS (.segment_dilated_sky from v9)
##      instead of a plain median on segim==0 pixels inside the stamp;
##    - GAMA PSFQUAL bit 64 when the median contributing segment size caps the
##      PSF postage stamp smaller than the box needed for target_flux_frac.
##
##  Source profuse_psf_v9.R BEFORE this file so the helpers exist (the
##  notebook build_psfs_v11.ipynb does this).
##
##  Method (per the GAMA DR4 BDInputsRv05 description and the paper's
##  Appendix A):
##
##    1. Star pool from ProFound segstats:
##         - edge fraction = 1 (no other-object / border / mask contact)
##         - edge excess < 1 (regular boundary geometry)
##         - axial ratio > 0.5
##         - not flagged as possibly spurious
##         - of these, joint cut in R50 (size) AND R50/R90 (concentration)
##           identifies the stellar locus (top fraction by concentration
##           inside the small-R50 group; default 4% as in GAMA gri stacks).
##         - bright > saturation_offset and faint < detection_offset of
##           the brightest star (5-sigma point source detection limit
##           and saturation are set per-band by KiDS/VIKING release notes;
##           we use offsets relative to the brightest clean star).
##
##    2. For each candidate, take a square cutout (default ~6*FWHM half-
##       width) and require:
##         - <10% of pixels belong to OTHER segments
##         - cutout fully inside the image
##         - sum(stamp) > 0
##       Then mask other segments, subtract a local sky from the segim==0
##       pixels of the stamp, and renormalise the cutout to magnitude 0
##       using the segstats flux of THIS star (sum_in_seg = 10^(-0.4*0)).
##
##    3. Per-star Moffat fit using ProFit::profitMakeModel and BFGS via
##       stats::optim(method="BFGS"), with:
##         - all 7 parameters free EXCEPT box (boxyness); fitted: xcen,
##           ycen, mag, FWHM, con, axrat, ang
##         - scale parameters fitted in log-space (log(FWHM) and
##           log(con-1) so that con > 1 is enforced)
##         - axial ratio fitted via logit so 0 < axrat <= 1
##         - normal likelihood: -0.5 * sum( ((data - model)/sigma)^2 )
##           with a constant per-star sigma estimated from non-segment
##           pixels in the stamp
##         - initial guesses from segstats (FWHM ~ 2*R50, con=2.5,
##           axrat from segstats axrat, ang from segstats ang)
##
##    4. Suitability after fit:
##         - |xcen_fit - cutout_centre| < 1 px in x AND y
##         - |mag_fit - 0| < 0.1
##         - reduced chi-square < 3
##         - none of FWHM, con, axrat, ang at the fit limits, except
##           axrat may equal its UPPER limit of 1 (perfectly round)
##         - iterative 2-sigma clip in log-space (where appropriate) on
##           FWHM, con, axrat, ang and background
##
##    5. Spatial selection from the suitable stars (relative to the
##       galaxy at the image centre by convention here):
##         - closest 2 from each of the 4 quadrants (8 total)
##         - if a quadrant has fewer than 2, fill with the closest
##           remaining suitable stars from any other quadrant
##         - if the total is < 8, use all of them
##         - if 0 suitable, use the DEFAULT PSF
##
##    6. Build the model PSF from the MEDIAN of the fitted parameters
##       (FWHM, con, axrat, ang).  Centre is forced to the centre of the
##       PSF image and the magnitude is forced to 0 so the image sums to 1.
##       The PSF size is adjusted so that the analytic Moffat contains
##       at least target_flux_frac (default 0.99) of its flux, capped at
##       the median segment size of the contributing stars (so we do
##       not extrapolate beyond the data).  The four corner pixels
##       outside the inscribed circle are set to 0 to avoid a strictly
##       rectangular PSF, then the image is renormalised.
##
##    7. Default PSF fallback: a perfectly round Moffat with
##       FWHM = `default_fwhm_arcsec / pixscale_arcsec`, con = 2.1
##       (the median r-band concentration measured by Casura+2022).
##
##  The orchestrator returns the same shape as v8/v9 PLUS extra slots so
##  callers can A/B compare against the GAMA-style PSF without rebuilding:
##
##    list(psf,                  # GAMA median-parameter Moffat (RECOMMENDED)
##         psf_default,          # the default Moffat (always built)
##         psf_per_star_models,  # NstarSel-list of individual fitted Moffats
##         psf_empirical,        # neighbour-masked recentred median stack (v9)
##         psf_apodized,         # empirical core + Moffat wing (v9)
##         output_kind,
##         star_table,           # full per-star fit + suitability + selection
##         decision,
##         segim, stars,
##         sky_pedestal_subtracted,
##         psf_quality_flag)     # bitwise as in GAMA DR4
## ============================================================================

suppressPackageStartupMessages({
  library(Rfits); library(ProFound); library(ProFit)
  library(magicaxis)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

## ProFound segmap ID at rounded pixel position (target anchor for exclusions)
.segim_id_at_xy <- function(segim, xy) {
  if (is.null(segim)) return(NA_integer_)
  if (length(xy) < 2L || !all(is.finite(xy))) return(NA_integer_)
  nr <- nrow(segim); nc <- ncol(segim)
  ix <- max(1L, min(nr, as.integer(round(xy[[1]]))))
  iy <- max(1L, min(nc, as.integer(round(xy[[2]]))))
  as.integer(segim[ix, iy])
}

## ----------------------------------------------------------------------------
##  Bring in v9 helpers when this file is sourced AFTER profuse_psf_v9.R
##  (which the build_psfs_v11 notebook does).  We use v9's Lanczos shift
##  and segment-dilated empirical PSF builder for the diagnostic sidecar
##  outputs psf_empirical and psf_apodized.
## ----------------------------------------------------------------------------
.has_v9 <- function() exists("build_empirical_psf_v9", inherits = TRUE) &&
                     exists(".lanczos_shift", inherits = TRUE) &&
                     exists("fit_moffat_to_stack", inherits = TRUE) &&
                     exists("build_apodized_psf", inherits = TRUE)

## ----------------------------------------------------------------------------
##  Render a Moffat PSF on a `dim`-pixel grid centred at the image centre,
##  with magnitude = 0 (so total flux = 10^0 = 1) and the requested shape.
##  Uses ProFit::profitMakeModel with no convolution (psf = NULL).
## ----------------------------------------------------------------------------
.render_moffat <- function(fwhm, con, axrat = 1, ang = 0, dim = c(75L, 75L),
                           magzero = 0) {
  dim <- as.integer(dim); if (any(dim %% 2L == 0L)) dim <- dim + 1L
  cx <- (dim[1] + 1) / 2; cy <- (dim[2] + 1) / 2
  ml <- list(moffat = list(
    xcen  = cx, ycen = cy,
    mag   = 0,
    fwhm  = fwhm,
    con   = con,
    axrat = axrat,
    ang   = ang,
    box   = 0
  ))
  out <- ProFit::profitMakeModel(modellist = ml, magzero = magzero,
                                 dim = dim, psf = NULL)
  z <- if (is.list(out) && !is.null(out$z)) out$z else out
  s <- sum(z, na.rm = TRUE)
  if (is.finite(s) && s > 0) z <- z / s
  attr(z, "fwhm")  <- fwhm
  attr(z, "beta")  <- con
  attr(z, "axrat") <- axrat
  attr(z, "ang")   <- ang
  z
}

## ----------------------------------------------------------------------------
##  Smallest odd dim that contains target_flux_frac of a circular Moffat's
##  flux on its analytic radial profile, capped at max_dim.
##
##  The Moffat enclosed flux at radius r (circular) is
##    F(<r)/F_inf = 1 - (1 + (r/alpha)^2)^(1 - beta)
##  where alpha = (FWHM/2) / sqrt(2^(1/beta) - 1).  Solve for r at frac.
## ----------------------------------------------------------------------------
.moffat_size_for_flux <- function(fwhm, con, target_flux_frac = 0.99,
                                  max_dim = 151L, min_dim = 21L) {
  if (!is.finite(fwhm) || !is.finite(con) || con <= 1)
    return(as.integer(max_dim))
  alpha <- (fwhm / 2) / sqrt(2^(1 / con) - 1)
  r_frac <- alpha * sqrt((1 - target_flux_frac)^(1 / (1 - con)) - 1)
  ## image needs to span +/- r_frac; round up to odd integer
  d <- 2 * ceiling(r_frac) + 1L
  d <- max(min_dim, min(max_dim, d))
  if (d %% 2L == 0L) d <- d + 1L
  as.integer(d)
}

## ----------------------------------------------------------------------------
##  Set the four corner pixels (outside the inscribed circle) to 0 so the
##  PSF image is not strictly rectangular.  Renormalise to sum=1.
## ----------------------------------------------------------------------------
.zero_corners <- function(psf) {
  if (is.null(psf)) return(psf)
  nr <- nrow(psf); nc <- ncol(psf)
  cx <- (nr + 1) / 2; cy <- (nc + 1) / 2
  rr <- sqrt(outer((seq_len(nr) - cx)^2, (seq_len(nc) - cy)^2, FUN = "+"))
  inscribed <- min(nr, nc) / 2
  out <- psf
  out[rr > inscribed] <- 0
  s <- sum(out, na.rm = TRUE)
  if (is.finite(s) && s > 0) out <- out / s
  attributes(out) <- attributes(psf)
  out
}

## ----------------------------------------------------------------------------
##  Iterative 2-sigma clip on a numeric vector.  In_log=TRUE clips on
##  the log of x (for scale parameters).  Returns a logical mask of the
##  KEPT entries (same length as x).
## ----------------------------------------------------------------------------
.iter_2sig_clip <- function(x, in_log = FALSE, sigma = 2, max_iter = 5) {
  ok <- is.finite(x)
  if (in_log) ok <- ok & x > 0
  if (sum(ok) < 3L) return(ok)
  v <- if (in_log) log(x[ok]) else x[ok]
  keep <- rep(TRUE, length(v))
  for (it in seq_len(max_iter)) {
    m <- stats::median(v[keep]); s <- stats::mad(v[keep]) * 1.4826
    if (!is.finite(s) || s == 0) break
    new_keep <- abs(v - m) < sigma * s
    if (identical(new_keep, keep)) break
    keep <- new_keep
    if (sum(keep) < 3L) break
  }
  out <- ok
  out[which(ok)] <- keep
  out
}

## ----------------------------------------------------------------------------
##  GAMA-style star pool selection from ProFound segstats.
##
##    pf            : ProFound output
##    n_max         : max stars to fit
##    edge_excess_max, axrat_min : hygiene cuts
##    star_frac     : fraction of clean candidates classed as star
##                    candidates by the joint R50 / R50R90 cut (4% in GAMA
##                    for gri stacks at skycut=2; tune for other depths)
##    sat_off, faint_off : magnitude offsets relative to the brightest
##                    accepted star (proxy for saturation / 5-sigma
##                    detection limit when KiDS/VIKING release-note values
##                    are not available locally)
##    exclude_segids   : drop these ProFound segIDs from the pool (e.g. the
##                    science target at the cutout centre is often *not* a star)
##
##  Returns a data.frame with one row per candidate (segID, segstats
##  values, R50R90).
## ----------------------------------------------------------------------------
select_clean_stars_v10 <- function(pf,
                                   n_max         = 50L,
                                   edge_excess_max = 1.0,
                                   axrat_min     = 0.5,
                                   star_frac     = 0.04,
                                   sat_off       = 0.5,
                                   faint_off     = 5.0,
                                   exclude_segids = NULL,
                                   verbose       = TRUE) {
  ss <- pf$segstats
  if (is.null(ss) || nrow(ss) == 0L) return(NULL)
  has_col <- function(nm) nm %in% names(ss)

  ## --- hygiene cuts -------------------------------------------------------
  no_touch <- ((ss$Nobject %||% 0) == 0) &
              ((ss$Nborder %||% 0) == 0) &
              ((ss$Nmask   %||% 0) == 0)
  edge_excess_ok <- if (has_col("edge_excess")) {
    is.finite(ss$edge_excess) & ss$edge_excess < edge_excess_max
  } else TRUE
  axrat_ok <- if (has_col("axrat")) {
    is.finite(ss$axrat) & ss$axrat > axrat_min
  } else TRUE
  ## ProFound attaches a sparse-source / quality flag column; clean = 0
  ## (must compare element-wise; isTRUE(...) on a vector is wrong here).
  flag_ok <- if (has_col("flag")) {
    is.na(ss$flag) | ss$flag == 0L
  } else TRUE
  if (!all(c("mag", "R50") %in% names(ss))) return(NULL)
  base <- no_touch & edge_excess_ok & axrat_ok & flag_ok &
          is.finite(ss$mag) & is.finite(ss$R50) & ss$R50 > 0
  if (has_col("R90")) base <- base & is.finite(ss$R90) & ss$R90 > ss$R50
  base[is.na(base)] <- FALSE
  if (!any(base)) {
    if (verbose) message("[psf-v11] no clean candidates after hygiene cuts.")
    return(NULL)
  }

  ## --- joint R50 / R50R90 cut: identify stellar locus ---------------------
  R50R90 <- if (has_col("R90")) ss$R50 / ss$R90 else NA_real_
  bright <- base
  ## use only the bright half for the locus to avoid the extended-source tail
  med_mag <- stats::median(ss$mag[base], na.rm = TRUE)
  if (sum(base & ss$mag < med_mag) >= 5) {
    bright <- base & ss$mag < med_mag
  }
  ## "stars are small AND highly concentrated":
  ## small = lower star_frac (e.g. 4%) of bright in R50;
  ## concentrated = lower star_frac of bright in R50/R90.
  R50_thr <- if (sum(bright) >= 1) {
    stats::quantile(ss$R50[bright], probs = star_frac, na.rm = TRUE)
  } else NA_real_
  conc_thr <- if (has_col("R90") && sum(bright) >= 1) {
    stats::quantile(R50R90[bright], probs = star_frac, na.rm = TRUE)
  } else NA_real_
  star_loc <- base &
              (if (is.finite(R50_thr))  ss$R50  <= R50_thr  else TRUE) &
              (if (is.finite(conc_thr)) R50R90 <= conc_thr else TRUE)
  star_loc[is.na(star_loc)] <- FALSE
  if (!any(star_loc)) {
    if (verbose) message("[psf-v11] joint R50/R50R90 cut left 0 candidates; ",
                         "loosening to base hygiene set.")
    star_loc <- base
  }

  ## --- saturation / faintness cuts (relative to brightest accepted) -------
  brt <- min(ss$mag[star_loc], na.rm = TRUE)
  star_loc <- star_loc & ss$mag >= brt + sat_off & ss$mag <= brt + faint_off
  if (!any(star_loc)) {
    if (verbose) message("[psf-v11] saturation/faintness cut emptied the pool.")
    return(NULL)
  }

  idx <- which(star_loc)
  if (!is.null(exclude_segids) && length(idx) > 0L) {
    ex <- unique(as.integer(exclude_segids[!is.na(exclude_segids)]))
    ex <- ex[is.finite(ex) & ex > 0L]
    if (length(ex))
      idx <- idx[!(ss$segID[idx] %in% ex)]
  }
  if (!length(idx)) {
    if (verbose)
      message("[psf-v11] star pool empty after segID exclusion ",
              "(target / user exclude_segids list).")
    return(NULL)
  }
  if (length(idx) > n_max) idx <- idx[order(ss$mag[idx])][seq_len(n_max)]
  out <- ss[idx, , drop = FALSE]
  out$R50R90 <- R50R90[idx]
  if (verbose) message(sprintf(
    "[psf-v11] star pool: %d candidates (R50<=%.2fpx, R50/R90<=%.3f, mag in [%.2f, %.2f])",
    nrow(out),
    R50_thr  %||% NA_real_,
    conc_thr %||% NA_real_,
    brt + sat_off, brt + faint_off))
  out
}

## ----------------------------------------------------------------------------
##  Per-star Moffat fit with ProFit + BFGS in log-space, GAMA-style.
##
##    image      : full image
##    segim      : full segmap
##    sky_rms    : full noise estimate (per-pixel sigma); NULL -> estimated
##                 from non-segment pixels in each stamp
##    seg_row    : one row of the candidate data.frame (segstats subset)
##    half_size  : half-width of the per-star cutout in pixels (the stamp
##                 will be (2*half_size + 1) on a side)
##    other_seg_pad : pad other segments by this many pixels before
##                 masking (covers undetected stellar haloes)
##    bfgs_iters : BFGS iterations
##    par_init   : optional initial-parameter override
##
##  Returns a list with:
##    ok, segID, stamp_size, segstats_mag, fit_par (named numeric),
##    fit_par_log (raw optim() solution), reduced_chisq, bg, bg_init,
##    fit_limits (logical for at-limit on each param), opt_convergence,
##    nllh, ndof
## ----------------------------------------------------------------------------
.fit_one_star_gama_v11 <- function(image, segim, seg_row,
                               half_size = 12L,
                               other_seg_pad = 1L,
                               sigma_floor = NULL,
                               bfgs_iters = 200L,
                               magzero = 0,
                               par_init = NULL,
                               other_pix_frac_max = 0.10,
                               positive_sum_check = TRUE,
                               return_stamp = FALSE,
                               verbose = FALSE) {
  nr <- nrow(image); nc <- ncol(image)
  segID <- as.integer(seg_row$segID)
  xc <- as.numeric(seg_row$xcen); yc <- as.numeric(seg_row$ycen)
  if (!is.finite(xc) || !is.finite(yc))
    return(list(ok = FALSE, segID = segID, reason = "bad_centroid"))
  ## v11: IRAF-style sub-pixel peak (v9) — reduces centroid bias before the stamp.
  if (exists(".refine_centroid", inherits = TRUE)) {
    xy <- tryCatch(.refine_centroid(image, xc, yc, win = 3L),
                   error = function(e) NULL)
    if (!is.null(xy) && length(xy) == 2L && all(is.finite(xy))) {
      xc <- xy[[1L]]; yc <- xy[[2L]]
    }
  }

  ## ------ Cutout ----------------------------------------------------------
  xri <- as.integer(round(xc)); yri <- as.integer(round(yc))
  hx <- as.integer(half_size); hy <- as.integer(half_size)
  x_lo <- xri - hx; x_hi <- xri + hx
  y_lo <- yri - hy; y_hi <- yri + hy
  if (x_lo < 1L || y_lo < 1L || x_hi > nr || y_hi > nc)
    return(list(ok = FALSE, segID = segID, reason = "near_image_edge"))

  data_stamp <- image[x_lo:x_hi, y_lo:y_hi, drop = FALSE]
  seg_stamp  <- if (is.null(segim)) matrix(0L, nrow(data_stamp), ncol(data_stamp))
                else segim[x_lo:x_hi, y_lo:y_hi, drop = FALSE]
  if (any(!is.finite(data_stamp))) data_stamp[!is.finite(data_stamp)] <- 0

  ## ------ Other-segment fraction check -----------------------------------
  other <- (seg_stamp != 0L) & (seg_stamp != segID)
  npx <- length(seg_stamp)
  if (sum(other) / npx > other_pix_frac_max)
    return(list(ok = FALSE, segID = segID, reason = "too_many_other_seg_pixels"))

  ## ------ Local sky: GAMA uses robust sky on the cutout; v11 prefers the
  ##        segment-dilated ring (v9) when profuse_psf_v9.R was sourced.
  sky_pix <- as.numeric(data_stamp[seg_stamp == 0L])
  bg_init <- 0
  used_seg_sky <- FALSE
  sky_ring_rms <- NA_real_
  if (exists(".segment_dilated_sky", inherits = TRUE) && !is.null(segim)) {
    fwhm0_est <- max(1.5, 2 * (seg_row$R50 %||% 1))
    fb_in  <- 3 * fwhm0_est
    fb_out <- max(fb_in + 2, 6 * fwhm0_est)
    sk <- tryCatch(
      .segment_dilated_sky(data_stamp, seg_stamp, segID,
                           inner_d = 2L, outer_d = 4L,
                           other_pad = other_seg_pad,
                           min_sky_pixels = 30L,
                           fallback_in = fb_in, fallback_out = fb_out),
      error = function(e) NULL)
    if (!is.null(sk) && is.finite(sk$sky)) {
      bg_init <- sk$sky
      used_seg_sky <- TRUE
      if (is.finite(sk$rms) && sk$rms > 0) sky_ring_rms <- sk$rms
    }
  }
  if (!used_seg_sky) {
    bg_init <- if (length(sky_pix) >= 8L) stats::median(sky_pix, na.rm = TRUE) else 0
  }
  if (!is.finite(bg_init)) bg_init <- 0
  data_stamp <- data_stamp - bg_init
  if (is.finite(sky_ring_rms) && sky_ring_rms > 0) {
    sigma_pix <- sky_ring_rms
  } else {
    sky_pix2 <- as.numeric(data_stamp[seg_stamp == 0L])
    sigma_pix <- if (length(sky_pix2) >= 8L) {
      stats::mad(sky_pix2, na.rm = TRUE) * 1.4826
    } else {
      stats::mad(as.numeric(data_stamp), na.rm = TRUE) * 1.4826
    }
  }
  if (!is.finite(sigma_pix) || sigma_pix <= 0) sigma_pix <- 1
  if (!is.null(sigma_floor)) sigma_pix <- max(sigma_pix, sigma_floor)

  ## ------ Renormalise stamp to magnitude 0 using segstats flux ------------
  star_pix <- (seg_stamp == segID)
  in_seg_flux <- sum(data_stamp[star_pix], na.rm = TRUE)
  if (positive_sum_check && (!is.finite(in_seg_flux) || in_seg_flux <= 0))
    return(list(ok = FALSE, segID = segID, reason = "non_positive_in_seg_flux"))
  if (!is.finite(in_seg_flux) || in_seg_flux <= 0) in_seg_flux <- 1
  norm_stamp <- data_stamp / in_seg_flux
  norm_sigma <- sigma_pix / in_seg_flux

  ## ------ Other-segment mask (pad by other_seg_pad) ----------------------
  fit_mask <- !other
  if (other_seg_pad > 0L) {
    pad <- matrix(FALSE, nrow(other), ncol(other))
    nb <- nrow(other); nc2 <- ncol(other)
    for (k in seq_len(other_seg_pad)) {
      pad <- other
      if (nb > 1L) {
        pad[1:(nb - 1L), ] <- pad[1:(nb - 1L), ] | other[2:nb,         ]
        pad[2:nb,         ] <- pad[2:nb,         ] | other[1:(nb - 1L), ]
      }
      if (nc2 > 1L) {
        pad[, 1:(nc2 - 1L)] <- pad[, 1:(nc2 - 1L)] | other[, 2:nc2        ]
        pad[, 2:nc2        ] <- pad[, 2:nc2        ] | other[, 1:(nc2 - 1L)]
      }
      other <- pad
    }
    fit_mask <- !other
  }
  ## also drop non-finite pixels
  fit_mask <- fit_mask & is.finite(norm_stamp)
  if (sum(fit_mask) < 25L)
    return(list(ok = FALSE, segID = segID, reason = "too_few_unmasked_pixels"))

  ## ------ Initial guesses from segstats ----------------------------------
  fwhm0  <- max(1.5, 2 * (seg_row$R50 %||% 1))
  con0   <- 2.5
  axrat0 <- max(min(seg_row$axrat %||% 0.95, 0.99), 0.2)
  ang0   <- (seg_row$ang %||% 0) %% 180
  ## centre of stamp in stamp coords
  ccx <- (nrow(norm_stamp) + 1) / 2
  ccy <- (ncol(norm_stamp) + 1) / 2
  par0 <- if (!is.null(par_init)) par_init else c(
    log_fwhm   = log(fwhm0),
    log_conm1  = log(max(con0 - 1, 0.1)),
    mag        = 0,
    axrat_logit = log(axrat0 / (1 - min(axrat0, 0.999))),
    ang        = ang0,
    xcen       = ccx,
    ycen       = ccy,
    bg         = 0
  )
  ## fit limits in PARAMETER space (for the at-limit check later)
  fit_lims <- list(
    fwhm  = c(0.5, 4 * fwhm0 + 5),
    con   = c(1.05, 20),
    axrat = c(0.1, 1.0),
    ang   = c(-180, 180),
    xcen  = c(ccx - 3, ccx + 3),
    ycen  = c(ccy - 3, ccy + 3),
    mag   = c(-1, 1)
  )

  unpack_par <- function(p) {
    list(
      fwhm  = exp(p[1]),
      con   = 1 + exp(p[2]),
      mag   = p[3],
      axrat = 1 / (1 + exp(-p[4])),
      ang   = p[5],
      xcen  = p[6],
      ycen  = p[7],
      bg    = p[8]
    )
  }

  build_model <- function(par) {
    pp <- unpack_par(par)
    if (!is.finite(pp$fwhm) || !is.finite(pp$con) ||
        pp$fwhm <= 0 || pp$con <= 1) return(NULL)
    ml <- list(moffat = list(
      xcen  = pp$xcen, ycen = pp$ycen,
      mag   = pp$mag,
      fwhm  = pp$fwhm,
      con   = pp$con,
      axrat = pp$axrat,
      ang   = pp$ang,
      box   = 0
    ))
    out <- tryCatch(
      ProFit::profitMakeModel(modellist = ml, magzero = magzero,
                              dim = dim(norm_stamp), psf = NULL),
      error = function(e) NULL)
    if (is.null(out)) return(NULL)
    z <- if (is.list(out) && !is.null(out$z)) out$z else out
    z + pp$bg
  }

  ## quadratic-soft-wall penalty so BFGS gets a smooth gradient near the
  ## boundary instead of a flat 1e6 cliff (which can confuse the line search)
  .soft_wall <- function(val, lo, hi, scale = 1e3) {
    pen <- 0
    if (val < lo) pen <- pen + scale * (lo - val)^2
    if (val > hi) pen <- pen + scale * (val - hi)^2
    pen
  }
  nllh <- function(par) {
    pp <- unpack_par(par)
    pen <- 0
    pen <- pen + .soft_wall(pp$fwhm, fit_lims$fwhm[1], fit_lims$fwhm[2])
    pen <- pen + .soft_wall(pp$con,  fit_lims$con[1],  fit_lims$con[2])
    pen <- pen + .soft_wall(pp$xcen, fit_lims$xcen[1], fit_lims$xcen[2])
    pen <- pen + .soft_wall(pp$ycen, fit_lims$ycen[1], fit_lims$ycen[2])
    pen <- pen + .soft_wall(pp$mag,  fit_lims$mag[1],  fit_lims$mag[2])
    pen <- pen + .soft_wall(pp$axrat, fit_lims$axrat[1], fit_lims$axrat[2])
    ## angle is circular but we still want to discourage extreme excursions
    ang_norm <- ((pp$ang + 180) %% 360) - 180
    pen <- pen + .soft_wall(ang_norm, fit_lims$ang[1], fit_lims$ang[2])
    if (pen > 1e8) return(pen)  # well outside; skip the model evaluation
    m <- build_model(par)
    if (is.null(m)) return(pen + 1e10)
    r <- (norm_stamp - m) / norm_sigma
    val <- 0.5 * sum(r[fit_mask]^2)
    if (!is.finite(val)) return(pen + 1e10)
    val + pen
  }

  fit <- tryCatch(
    stats::optim(par0, nllh, method = "BFGS",
                 control = list(maxit = bfgs_iters, reltol = 1e-7)),
    error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$value))
    return(list(ok = FALSE, segID = segID, reason = "bfgs_failed"))

  pp <- unpack_par(fit$par)
  m  <- build_model(fit$par)
  if (is.null(m))
    return(list(ok = FALSE, segID = segID, reason = "model_eval_failed"))
  resid <- (norm_stamp - m) / norm_sigma
  ndof <- max(1L, sum(fit_mask) - length(fit$par))
  rchisq <- sum(resid[fit_mask]^2) / ndof

  ## at-limit check (1% relative tolerance for scale params, 0.05 for axrat)
  at_lim <- list(
    fwhm  = pp$fwhm  <= fit_lims$fwhm[1]  * 1.01 ||
            pp$fwhm  >= fit_lims$fwhm[2]  * 0.99,
    con   = pp$con   <= fit_lims$con[1]   * 1.01 ||
            pp$con   >= fit_lims$con[2]   * 0.99,
    axrat = pp$axrat <= fit_lims$axrat[1] * 1.01,   # axrat upper limit (1) is OK
    ang   = abs(((pp$ang + 180) %% 360 - 180) - fit_lims$ang[1]) < 0.5 ||
            abs(((pp$ang + 180) %% 360 - 180) - fit_lims$ang[2]) < 0.5
  )

  list(
    ok = TRUE, segID = segID,
    stamp_size = nrow(norm_stamp),
    segstats_mag = seg_row$mag,
    segstats_R50 = seg_row$R50,
    segstats_R90 = if ("R90" %in% names(seg_row)) seg_row$R90 else NA_real_,
    segstats_xcen = xc, segstats_ycen = yc,
    fit_par = pp,
    fit_par_log = fit$par,
    bg_init = bg_init,
    sigma_pix = sigma_pix,
    nllh = fit$value, ndof = ndof, reduced_chisq = rchisq,
    at_lim = at_lim,
    convergence = fit$convergence,
    stamp = if (return_stamp) norm_stamp else NULL,
    fit_mask = if (return_stamp) fit_mask else NULL,
    cutout_anchor = c(x_lo = x_lo, y_lo = y_lo, x_hi = x_hi, y_hi = y_hi)
  )
}

## ----------------------------------------------------------------------------
##  v12: per-star Moffat fit using L-BFGS-B with hard bounds.
##
##  Drop-in replacement for .fit_one_star_gama_v11.  Identical inputs and
##  outputs.  The only changes versus v11:
##    1. `optim(method="L-BFGS-B", lower=, upper=)` instead of BFGS + soft
##       walls.  The optimiser cannot leave the box.
##    2. The fit_lims for `mag` are tightened from c(-1, 1) to c(-0.3, 0.3)
##       since the stamp is renormalised so the in-segment flux is exactly
##       1.0 (mag=0 by construction).  A real star can drift by < 0.3 mag
##       from the segstats mag; anything beyond that is the optimiser
##       running away.
##    3. `bg` is now bounded to +/- 5 * sigma_pix (was unbounded).
##    4. Initial bg is the segment-dilated sky pre-subtraction value.
##
## ----------------------------------------------------------------------------
.fit_one_star_gama_v12 <- function(image, segim, seg_row,
                               half_size = 12L,
                               other_seg_pad = 1L,
                               sigma_floor = NULL,
                               bfgs_iters = 200L,
                               magzero = 0,
                               par_init = NULL,
                               other_pix_frac_max = 0.10,
                               positive_sum_check = TRUE,
                               return_stamp = FALSE,
                               verbose = FALSE) {
  nr <- nrow(image); nc <- ncol(image)
  segID <- as.integer(seg_row$segID)
  xc <- as.numeric(seg_row$xcen); yc <- as.numeric(seg_row$ycen)
  if (!is.finite(xc) || !is.finite(yc))
    return(list(ok = FALSE, segID = segID, reason = "bad_centroid"))
  if (exists(".refine_centroid", inherits = TRUE)) {
    xy <- tryCatch(.refine_centroid(image, xc, yc, win = 3L),
                   error = function(e) NULL)
    if (!is.null(xy) && length(xy) == 2L && all(is.finite(xy))) {
      xc <- xy[[1L]]; yc <- xy[[2L]]
    }
  }

  xri <- as.integer(round(xc)); yri <- as.integer(round(yc))
  hx <- as.integer(half_size); hy <- as.integer(half_size)
  x_lo <- xri - hx; x_hi <- xri + hx
  y_lo <- yri - hy; y_hi <- yri + hy
  if (x_lo < 1L || y_lo < 1L || x_hi > nr || y_hi > nc)
    return(list(ok = FALSE, segID = segID, reason = "near_image_edge"))

  data_stamp <- image[x_lo:x_hi, y_lo:y_hi, drop = FALSE]
  seg_stamp  <- if (is.null(segim)) matrix(0L, nrow(data_stamp), ncol(data_stamp))
                else segim[x_lo:x_hi, y_lo:y_hi, drop = FALSE]
  if (any(!is.finite(data_stamp))) data_stamp[!is.finite(data_stamp)] <- 0

  other <- (seg_stamp != 0L) & (seg_stamp != segID)
  npx <- length(seg_stamp)
  if (sum(other) / npx > other_pix_frac_max)
    return(list(ok = FALSE, segID = segID, reason = "too_many_other_seg_pixels"))

  ## ------ Local sky from segment-dilated ring (or robust median fallback) -
  sky_pix <- as.numeric(data_stamp[seg_stamp == 0L])
  bg_init <- 0
  used_seg_sky <- FALSE
  sky_ring_rms <- NA_real_
  if (exists(".segment_dilated_sky", inherits = TRUE) && !is.null(segim)) {
    fwhm0_est <- max(1.5, 2 * (seg_row$R50 %||% 1))
    fb_in  <- 3 * fwhm0_est
    fb_out <- max(fb_in + 2, 6 * fwhm0_est)
    sk <- tryCatch(
      .segment_dilated_sky(data_stamp, seg_stamp, segID,
                           inner_d = 2L, outer_d = 4L,
                           other_pad = other_seg_pad,
                           min_sky_pixels = 30L,
                           fallback_in = fb_in, fallback_out = fb_out),
      error = function(e) NULL)
    if (!is.null(sk) && is.finite(sk$sky)) {
      bg_init <- sk$sky
      used_seg_sky <- TRUE
      if (is.finite(sk$rms) && sk$rms > 0) sky_ring_rms <- sk$rms
    }
  }
  if (!used_seg_sky) {
    bg_init <- if (length(sky_pix) >= 8L) stats::median(sky_pix, na.rm = TRUE) else 0
  }
  if (!is.finite(bg_init)) bg_init <- 0
  data_stamp <- data_stamp - bg_init
  if (is.finite(sky_ring_rms) && sky_ring_rms > 0) {
    sigma_pix <- sky_ring_rms
  } else {
    sky_pix2 <- as.numeric(data_stamp[seg_stamp == 0L])
    sigma_pix <- if (length(sky_pix2) >= 8L) {
      stats::mad(sky_pix2, na.rm = TRUE) * 1.4826
    } else {
      stats::mad(as.numeric(data_stamp), na.rm = TRUE) * 1.4826
    }
  }
  if (!is.finite(sigma_pix) || sigma_pix <= 0) sigma_pix <- 1
  if (!is.null(sigma_floor)) sigma_pix <- max(sigma_pix, sigma_floor)

  ## ------ Renormalise to magnitude 0 ---------------------------------------
  star_pix <- (seg_stamp == segID)
  in_seg_flux <- sum(data_stamp[star_pix], na.rm = TRUE)
  if (positive_sum_check && (!is.finite(in_seg_flux) || in_seg_flux <= 0))
    return(list(ok = FALSE, segID = segID, reason = "non_positive_in_seg_flux"))
  if (!is.finite(in_seg_flux) || in_seg_flux <= 0) in_seg_flux <- 1
  norm_stamp <- data_stamp / in_seg_flux
  norm_sigma <- sigma_pix / in_seg_flux

  fit_mask <- !other
  if (other_seg_pad > 0L) {
    nb <- nrow(other); nc2 <- ncol(other)
    for (k in seq_len(other_seg_pad)) {
      pad <- other
      if (nb > 1L) {
        pad[1:(nb - 1L), ] <- pad[1:(nb - 1L), ] | other[2:nb,         ]
        pad[2:nb,         ] <- pad[2:nb,         ] | other[1:(nb - 1L), ]
      }
      if (nc2 > 1L) {
        pad[, 1:(nc2 - 1L)] <- pad[, 1:(nc2 - 1L)] | other[, 2:nc2        ]
        pad[, 2:nc2        ] <- pad[, 2:nc2        ] | other[, 1:(nc2 - 1L)]
      }
      other <- pad
    }
    fit_mask <- !other
  }
  fit_mask <- fit_mask & is.finite(norm_stamp)
  if (sum(fit_mask) < 25L)
    return(list(ok = FALSE, segID = segID, reason = "too_few_unmasked_pixels"))

  ## ------ Initial guesses --------------------------------------------------
  fwhm0  <- max(1.5, 2 * (seg_row$R50 %||% 1))
  con0   <- 2.5
  axrat0 <- max(min(seg_row$axrat %||% 0.95, 0.99), 0.2)
  ang0   <- (seg_row$ang %||% 0) %% 180
  ccx <- (nrow(norm_stamp) + 1) / 2
  ccy <- (ncol(norm_stamp) + 1) / 2
  par0 <- if (!is.null(par_init)) par_init else c(
    log_fwhm   = log(fwhm0),
    log_conm1  = log(max(con0 - 1, 0.1)),
    mag        = 0,
    axrat_logit = log(axrat0 / (1 - min(axrat0, 0.999))),
    ang        = ang0,
    xcen       = ccx,
    ycen       = ccy,
    bg         = 0
  )

  fit_lims <- list(
    fwhm  = c(0.5, 4 * fwhm0 + 5),
    con   = c(1.05, 20),
    axrat = c(0.1, 1.0),
    ang   = c(-180, 180),
    xcen  = c(ccx - 3, ccx + 3),
    ycen  = c(ccy - 3, ccy + 3),
    mag   = c(-0.3, 0.3)            # v12: tightened from c(-1, 1)
  )
  ## bg gets +/- 5 * norm_sigma room around its initial value (which is 0
  ## post-subtraction, so this is +/- 5 sigma in normalised flux units).
  bg_halfw <- max(5 * norm_sigma, 1e-6)

  ## v12: HARD bounds in TRANSFORMED-PARAMETER space (the space optim sees)
  lower <- c(
    log_fwhm   = log(fit_lims$fwhm[1]),
    log_conm1  = log(fit_lims$con[1] - 1),
    mag        = fit_lims$mag[1],
    axrat_logit = log(fit_lims$axrat[1] / (1 - fit_lims$axrat[1])),
    ang        = fit_lims$ang[1],
    xcen       = fit_lims$xcen[1],
    ycen       = fit_lims$ycen[1],
    bg         = -bg_halfw
  )
  upper <- c(
    log_fwhm   = log(fit_lims$fwhm[2]),
    log_conm1  = log(fit_lims$con[2] - 1),
    mag        = fit_lims$mag[2],
    axrat_logit = log(0.999 / (1 - 0.999)),     # approach axrat=1 smoothly
    ang        = fit_lims$ang[2],
    xcen       = fit_lims$xcen[2],
    ycen       = fit_lims$ycen[2],
    bg         =  bg_halfw
  )
  ## clamp par0 inside the box
  par0 <- pmin(pmax(par0, lower), upper)

  unpack_par <- function(p) {
    list(
      fwhm  = exp(p[1]),
      con   = 1 + exp(p[2]),
      mag   = p[3],
      axrat = 1 / (1 + exp(-p[4])),
      ang   = p[5],
      xcen  = p[6],
      ycen  = p[7],
      bg    = p[8]
    )
  }

  build_model <- function(par) {
    pp <- unpack_par(par)
    if (!is.finite(pp$fwhm) || !is.finite(pp$con) ||
        pp$fwhm <= 0 || pp$con <= 1) return(NULL)
    ml <- list(moffat = list(
      xcen  = pp$xcen, ycen = pp$ycen,
      mag   = pp$mag,
      fwhm  = pp$fwhm,
      con   = pp$con,
      axrat = pp$axrat,
      ang   = pp$ang,
      box   = 0
    ))
    out <- tryCatch(
      ProFit::profitMakeModel(modellist = ml, magzero = magzero,
                              dim = dim(norm_stamp), psf = NULL),
      error = function(e) NULL)
    if (is.null(out)) return(NULL)
    z <- if (is.list(out) && !is.null(out$z)) out$z else out
    z + pp$bg
  }

  nllh <- function(par) {
    if (any(!is.finite(par))) return(.Machine$double.xmax / 1e6)
    m <- build_model(par)
    if (is.null(m)) return(.Machine$double.xmax / 1e6)
    r <- (norm_stamp - m) / norm_sigma
    val <- 0.5 * sum(r[fit_mask]^2)
    if (!is.finite(val)) return(.Machine$double.xmax / 1e6)
    val
  }

  fit <- tryCatch(
    stats::optim(par0, nllh, method = "L-BFGS-B",
                 lower = lower, upper = upper,
                 control = list(maxit = bfgs_iters,
                                factr = 1e7,
                                pgtol = 1e-7)),
    error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$value))
    return(list(ok = FALSE, segID = segID, reason = "bfgs_failed"))

  pp <- unpack_par(fit$par)
  m  <- build_model(fit$par)
  if (is.null(m))
    return(list(ok = FALSE, segID = segID, reason = "model_eval_failed"))
  resid <- (norm_stamp - m) / norm_sigma
  ndof <- max(1L, sum(fit_mask) - length(fit$par))
  rchisq <- sum(resid[fit_mask]^2) / ndof

  ## at-limit checks: now meaningful because L-BFGS-B uses HARD bounds in
  ## the transformed-parameter space.  We test at-limit by checking whether
  ## the transformed parameter is within 1% of its lower/upper bound.
  rel_at_lim <- function(val, lo, hi, tol = 0.01) {
    abs(val - lo) < tol * max(1, abs(lo)) ||
      abs(val - hi) < tol * max(1, abs(hi))
  }
  at_lim <- list(
    fwhm  = rel_at_lim(fit$par[1], lower[1], upper[1]),
    con   = rel_at_lim(fit$par[2], lower[2], upper[2]),
    axrat = pp$axrat <= fit_lims$axrat[1] * 1.01,  # axrat upper (1) is OK
    ang   = abs(((pp$ang + 180) %% 360 - 180) - fit_lims$ang[1]) < 0.5 ||
            abs(((pp$ang + 180) %% 360 - 180) - fit_lims$ang[2]) < 0.5
  )

  list(
    ok = TRUE, segID = segID,
    stamp_size = nrow(norm_stamp),
    segstats_mag = seg_row$mag,
    segstats_R50 = seg_row$R50,
    segstats_R90 = if ("R90" %in% names(seg_row)) seg_row$R90 else NA_real_,
    segstats_xcen = xc, segstats_ycen = yc,
    fit_par = pp,
    fit_par_log = fit$par,
    bg_init = bg_init,
    sigma_pix = sigma_pix,
    nllh = fit$value, ndof = ndof, reduced_chisq = rchisq,
    at_lim = at_lim,
    convergence = fit$convergence,
    stamp = if (return_stamp) norm_stamp else NULL,
    fit_mask = if (return_stamp) fit_mask else NULL,
    cutout_anchor = c(x_lo = x_lo, y_lo = y_lo, x_hi = x_hi, y_hi = y_hi)
  )
}

## ----------------------------------------------------------------------------
##  GAMA-style suitability filter for the per-star fits.
##
##  Applies (in order):
##   - centre within +/- 1 px of stamp centre
##   - mag within +/- 0.1 of 0
##   - reduced chi-square < chisq_max
##   - none of FWHM, con, axrat, ang at the fit limits (axrat upper OK)
##   - iterative 2-sigma clip on log(FWHM), log(con-1), logit(axrat), ang, bg
##
##  Returns a logical vector "suitable" of length(fits).
## ----------------------------------------------------------------------------
.gama_suitability <- function(fits,
                              centre_tol_px = 1.0,
                              mag_tol_mag   = 0.1,
                              chisq_max     = 3,
                              do_2sig_clip  = TRUE,
                              verbose       = TRUE) {
  ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
  N <- length(fits)
  s <- rep(FALSE, N)
  reasons <- rep(NA_character_, N)
  for (i in which(ok)) {
    f <- fits[[i]]
    pp <- f$fit_par
    ccx <- (f$stamp_size + 1) / 2
    ccy <- ccx
    ok_centre <- abs(pp$xcen - ccx) <= centre_tol_px &&
                 abs(pp$ycen - ccy) <= centre_tol_px
    ok_mag    <- abs(pp$mag) <= mag_tol_mag
    ok_chisq  <- is.finite(f$reduced_chisq) && f$reduced_chisq < chisq_max
    ok_lims   <- !(isTRUE(f$at_lim$fwhm) || isTRUE(f$at_lim$con) ||
                   isTRUE(f$at_lim$axrat) || isTRUE(f$at_lim$ang))
    s[i] <- ok_centre && ok_mag && ok_chisq && ok_lims
    if (!s[i]) {
      reasons[i] <- paste0(
        if (!ok_centre) "centre " else "",
        if (!ok_mag)    "mag "    else "",
        if (!ok_chisq)  "chisq "  else "",
        if (!ok_lims)   "atlim "  else "")
    }
  }
  ## Iterative 2-sigma clip in log-space where appropriate.
  if (do_2sig_clip && sum(s) >= 5L) {
    fwhms <- vapply(fits, function(f) if (isTRUE(f$ok)) f$fit_par$fwhm  else NA_real_, numeric(1))
    conms <- vapply(fits, function(f) if (isTRUE(f$ok)) f$fit_par$con - 1 else NA_real_, numeric(1))
    axrs  <- vapply(fits, function(f) if (isTRUE(f$ok)) f$fit_par$axrat else NA_real_, numeric(1))
    angs  <- vapply(fits, function(f) if (isTRUE(f$ok)) ((f$fit_par$ang + 90) %% 180) - 90 else NA_real_, numeric(1))
    bgs   <- vapply(fits, function(f) if (isTRUE(f$ok)) f$fit_par$bg else NA_real_, numeric(1))
    keep <- s
    keep <- keep & .iter_2sig_clip(fwhms, in_log = TRUE)
    keep <- keep & .iter_2sig_clip(conms, in_log = TRUE)
    keep <- keep & .iter_2sig_clip(axrs,  in_log = TRUE)
    keep <- keep & .iter_2sig_clip(angs,  in_log = FALSE)
    keep <- keep & .iter_2sig_clip(bgs,   in_log = FALSE)
    n_clipped <- sum(s & !keep)
    s <- keep
    if (verbose && n_clipped > 0L)
      message(sprintf("[psf-v11] 2-sigma clip rejected %d additional star(s)",
                      n_clipped))
  }
  if (verbose) {
    n_failed <- sum(!ok)
    n_unsuit <- sum(ok & !s)
    message(sprintf(
      "[psf-v11] suitability: %d/%d stars suitable (failed_fit=%d, unsuitable=%d)",
      sum(s), N, n_failed, n_unsuit))
    if (n_unsuit > 0L) {
      bad <- which(ok & !s)
      bad <- bad[seq_len(min(8L, length(bad)))]
      for (j in bad)
        message(sprintf("  - segID=%s reason=%s",
                        fits[[j]]$segID, reasons[j] %||% "clip"))
    }
  }
  s
}

## BDInputsRv05: after basic gates, iterative 2-sigma rejection on **FWHM,
## concentration (con−1), axial ratio and angle only** — not sky/bg (Casura DR4).
.gama_suitability_dr4 <- function(fits,
                                  centre_tol_px = 1.0,
                                  mag_tol_mag   = 0.1,
                                  chisq_max     = 3,
                                  do_2sig_clip  = TRUE,
                                  verbose       = TRUE) {
  ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
  N <- length(fits)
  s <- rep(FALSE, N)
  reasons <- rep(NA_character_, N)
  for (i in which(ok)) {
    f <- fits[[i]]
    pp <- f$fit_par
    ccx <- (f$stamp_size + 1) / 2
    ccy <- ccx
    ok_centre <- abs(pp$xcen - ccx) <= centre_tol_px &&
                 abs(pp$ycen - ccy) <= centre_tol_px
    ok_mag    <- abs(pp$mag) <= mag_tol_mag
    ok_chisq  <- is.finite(f$reduced_chisq) && f$reduced_chisq < chisq_max
    ok_lims   <- !(isTRUE(f$at_lim$fwhm) || isTRUE(f$at_lim$con) ||
                   isTRUE(f$at_lim$axrat) || isTRUE(f$at_lim$ang))
    s[i] <- ok_centre && ok_mag && ok_chisq && ok_lims
    if (!s[i]) {
      reasons[i] <- paste0(
        if (!ok_centre) "centre " else "",
        if (!ok_mag)    "mag "    else "",
        if (!ok_chisq)  "chisq "  else "",
        if (!ok_lims)   "atlim "  else "")
    }
  }
  if (do_2sig_clip && sum(s) >= 5L) {
    fwhms <- vapply(fits, function(f)
      if (isTRUE(f$ok)) f$fit_par$fwhm else NA_real_, numeric(1))
    conms <- vapply(fits, function(f)
      if (isTRUE(f$ok)) f$fit_par$con - 1 else NA_real_, numeric(1))
    axrs  <- vapply(fits, function(f)
      if (isTRUE(f$ok)) f$fit_par$axrat else NA_real_, numeric(1))
    angs  <- vapply(fits, function(f) if (isTRUE(f$ok))
      ((f$fit_par$ang + 90) %% 180) - 90 else NA_real_, numeric(1))
    keep <- s
    keep <- keep & .iter_2sig_clip(fwhms, in_log = TRUE)
    keep <- keep & .iter_2sig_clip(conms, in_log = TRUE)
    keep <- keep & .iter_2sig_clip(axrs,  in_log = TRUE)
    keep <- keep & .iter_2sig_clip(angs,  in_log = FALSE)
    n_clipped <- sum(s & !keep)
    s <- keep
    if (verbose && n_clipped > 0L)
      message(sprintf("[psf-v13] BDInputs 2-sigma clip rejected %d star(s)",
                      n_clipped))
  }
  if (verbose) {
    message(sprintf(
      "[psf-v13] suitability (DR4): %d/%d stars suitable (failed=%d, excluded=%d)",
      sum(s), N, sum(!ok), sum(ok & !s)))
  }
  s
}

## ----------------------------------------------------------------------------
##  Quadrant-balanced selection of the closest stars to (gx, gy).
##
##  Returns up to 8 indices (relative to the suitable subset) chosen as:
##    - the 2 closest in each of the 4 quadrants (NE, NW, SE, SW)
##    - if a quadrant has < 2 suitable stars, the remaining slots are
##      filled by the closest unused suitable stars from any other
##      quadrant
##    - if there are < 8 in total, all are returned
##
##  Inputs:
##    star_xy : data.frame with xcen, ycen for each suitable star
##    gx, gy  : galaxy reference position (image coords)
## ----------------------------------------------------------------------------
.quadrant_select <- function(star_xy, gx, gy, n_per_quad = 2L, n_total = 8L) {
  if (nrow(star_xy) == 0L) return(integer(0))
  d <- sqrt((star_xy$xcen - gx)^2 + (star_xy$ycen - gy)^2)
  qx <- ifelse(star_xy$xcen >= gx, 1L, -1L)
  qy <- ifelse(star_xy$ycen >= gy, 1L, -1L)
  quads <- list()
  quads[["+x+y"]] <- which(qx == 1L  & qy == 1L)
  quads[["-x+y"]] <- which(qx == -1L & qy == 1L)
  quads[["+x-y"]] <- which(qx == 1L  & qy == -1L)
  quads[["-x-y"]] <- which(qx == -1L & qy == -1L)
  picked <- integer(0)
  for (q in names(quads)) {
    idx <- quads[[q]]
    if (length(idx) == 0L) next
    ord <- idx[order(d[idx])]
    picked <- c(picked, head(ord, n_per_quad))
  }
  ## fill remaining slots with the closest unused stars (if any quadrant
  ## was underfilled)
  if (length(picked) < n_total) {
    leftover <- setdiff(seq_len(nrow(star_xy)), picked)
    if (length(leftover) > 0L) {
      ord <- leftover[order(d[leftover])]
      need <- n_total - length(picked)
      picked <- c(picked, head(ord, need))
    }
  }
  if (length(picked) > n_total) picked <- picked[seq_len(n_total)]
  picked
}

## ----------------------------------------------------------------------------
##  Build the GAMA model PSF from the median of the fitted parameters.
##
##  Inputs:
##    selected_fits : list of suitable, selected per-star fits
##    target_flux_frac : 0.99 (Casura+2022)
##    max_dim : hard ceiling on the PSF size
##    median_segment_size : cap from the contributing stars' segments
##    pixscale_arcsec : if available, attached to attr
##
##  Returns a Moffat PSF with sum=1, corners zeroed, attrs (fwhm, beta,
##  axrat, ang, Nstar_used, fwhm_arcsec).
## ----------------------------------------------------------------------------
.build_median_param_psf <- function(selected_fits,
                                    target_flux_frac = 0.99,
                                    max_dim = 151L,
                                    min_dim = 21L,
                                    median_segment_size_px = NA_integer_,
                                    pixscale_arcsec = NA_real_) {
  if (length(selected_fits) == 0L) return(NULL)
  fwhms <- vapply(selected_fits, function(f) f$fit_par$fwhm,  numeric(1))
  cons  <- vapply(selected_fits, function(f) f$fit_par$con,   numeric(1))
  axrs  <- vapply(selected_fits, function(f) f$fit_par$axrat, numeric(1))
  angs  <- vapply(selected_fits, function(f) f$fit_par$ang,   numeric(1))

  fwhm_med <- stats::median(fwhms, na.rm = TRUE)
  con_med  <- stats::median(cons,  na.rm = TRUE)
  axr_med  <- stats::median(axrs,  na.rm = TRUE)
  ## angle: circular median via vector mean of cos(2*ang), sin(2*ang)
  ang_rad  <- 2 * angs * pi / 180
  ang_med  <- (atan2(stats::median(sin(ang_rad), na.rm = TRUE),
                     stats::median(cos(ang_rad), na.rm = TRUE))
              ) * 180 / pi / 2
  ang_med  <- ((ang_med + 180) %% 180)

  ## Size: contain target_flux_frac of the analytic Moffat, capped at the
  ## median segment size (so we never extrapolate beyond the data).
  d_flux <- .moffat_size_for_flux(fwhm_med, con_med, target_flux_frac,
                                  max_dim = max_dim, min_dim = min_dim)
  d_use  <- if (is.finite(median_segment_size_px) && median_segment_size_px > 0L)
              min(d_flux, as.integer(median_segment_size_px))
            else d_flux
  d_use  <- max(min_dim, d_use)
  if (d_use %% 2L == 0L) d_use <- d_use + 1L
  d_use  <- as.integer(d_use)

  psf <- .render_moffat(fwhm = fwhm_med, con = con_med, axrat = axr_med,
                        ang = ang_med, dim = c(d_use, d_use))
  psf <- .zero_corners(psf)

  attr(psf, "fwhm")  <- fwhm_med
  attr(psf, "beta")  <- con_med
  attr(psf, "axrat") <- axr_med
  attr(psf, "ang")   <- ang_med
  attr(psf, "Nstar_used") <- length(selected_fits)
  attr(psf, "kind")  <- "gama_median_param_moffat"
  attr(psf, "target_flux_frac") <- target_flux_frac
  attr(psf, "size_px") <- d_use
  if (is.finite(pixscale_arcsec)) {
    attr(psf, "pixscale_arcsec") <- pixscale_arcsec
    attr(psf, "fwhm_arcsec")     <- fwhm_med * pixscale_arcsec
  }
  psf
}

## ----------------------------------------------------------------------------
##  Build the default Moffat PSF (perfectly round) used when no suitable
##  stars are found.  Casura+2022 uses con = 2.1 (median r-band) and the
##  FWHM from the KiDS header.  We let the caller pass default_fwhm_px.
## ----------------------------------------------------------------------------
.build_default_psf <- function(default_fwhm_px, default_con = 2.1,
                               psf_dim = c(75L, 75L),
                               pixscale_arcsec = NA_real_) {
  d <- .moffat_size_for_flux(default_fwhm_px, default_con, 0.99,
                             max_dim = max(psf_dim), min_dim = 21L)
  psf <- .render_moffat(fwhm = default_fwhm_px, con = default_con,
                        axrat = 1, ang = 0, dim = c(d, d))
  psf <- .zero_corners(psf)
  attr(psf, "fwhm")  <- default_fwhm_px
  attr(psf, "beta")  <- default_con
  attr(psf, "axrat") <- 1
  attr(psf, "ang")   <- 0
  attr(psf, "kind")  <- "gama_default_round_moffat"
  attr(psf, "Nstar_used") <- 0L
  if (is.finite(pixscale_arcsec)) {
    attr(psf, "pixscale_arcsec") <- pixscale_arcsec
    attr(psf, "fwhm_arcsec")     <- default_fwhm_px * pixscale_arcsec
  }
  psf
}

## ----------------------------------------------------------------------------
##  build_psf_for_image_v11: v10 GAMA logic + v9 centroid/sky refinements.
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
##  v12: adaptive star-pool wrapper.  Calls v10 with progressively looser
##  cuts when the strict pool returns 0 candidates.  Reports which "tier"
##  was used so the caller can flag galaxies that needed the loose pool.
##
##  Tiers (in order):
##    1  strict          : star_frac=0.04, sat_off=0.5,  faint_off=5.0
##    2  loose-quantile  : star_frac=0.10, sat_off=0.5,  faint_off=5.0
##    3  loose-mag       : star_frac=0.10, sat_off=0.0,  faint_off=7.0
##    4  very-loose      : star_frac=0.25, sat_off=0.0,  faint_off=7.0
##    5  no-locus-cut    : star_frac=1.00, sat_off=0.0,  faint_off=8.0
##
##  Returns a list(cand, tier) where cand is the data.frame (or NULL) and
##  tier is the integer tier that produced it (or NA if all tiers failed).
##  exclude_segids is forwarded to select_clean_stars_v10 on every tier try.
## ----------------------------------------------------------------------------
select_clean_stars_v12 <- function(pf,
                                   n_max         = 50L,
                                   edge_excess_max = 1.0,
                                   axrat_min     = 0.5,
                                   star_frac     = 0.04,
                                   sat_off       = 0.5,
                                   faint_off     = 5.0,
                                   exclude_segids = NULL,
                                   verbose       = TRUE) {
  tiers <- list(
    list(star_frac = star_frac, sat_off = sat_off, faint_off = faint_off,
         label = "strict"),
    list(star_frac = max(star_frac, 0.10), sat_off = sat_off,
         faint_off = faint_off,             label = "loose-quantile"),
    list(star_frac = max(star_frac, 0.10), sat_off = 0.0,
         faint_off = max(faint_off, 7.0),   label = "loose-mag"),
    list(star_frac = max(star_frac, 0.25), sat_off = 0.0,
         faint_off = max(faint_off, 7.0),   label = "very-loose"),
    list(star_frac = 1.0,                  sat_off = 0.0,
         faint_off = max(faint_off, 8.0),   label = "no-locus-cut")
  )
  for (k in seq_along(tiers)) {
    t <- tiers[[k]]
    cand <- select_clean_stars_v10(
      pf,
      n_max         = n_max,
      edge_excess_max = edge_excess_max,
      axrat_min     = axrat_min,
      star_frac     = t$star_frac,
      sat_off       = t$sat_off,
      faint_off     = t$faint_off,
      exclude_segids = exclude_segids,
      verbose       = FALSE
    )
    if (!is.null(cand) && nrow(cand) > 0L) {
      if (verbose) message(sprintf(
        "[psf-v12] pool tier %d (%s): %d candidates (star_frac=%.2f sat_off=%.1f faint_off=%.1f)",
        k, t$label, nrow(cand), t$star_frac, t$sat_off, t$faint_off))
      return(list(cand = cand, tier = k, tier_label = t$label,
                  star_frac = t$star_frac, sat_off = t$sat_off,
                  faint_off = t$faint_off))
    }
    if (verbose) message(sprintf(
      "[psf-v12] pool tier %d (%s): empty, retrying...", k, t$label))
  }
  if (verbose) message("[psf-v12] all 5 pool tiers empty.")
  list(cand = NULL, tier = NA_integer_, tier_label = "none",
       star_frac = NA_real_, sat_off = NA_real_, faint_off = NA_real_)
}

build_psf_for_image_v11 <- function(image, magzero,
                                    pixscale_arcsec = NULL,
                                    psf_dim = c(75L, 75L),
                                    Nstar_pool = 50L,
                                    star_frac = 0.04,
                                    edge_excess_max = 1.0,
                                    star_pool_axrat_min = 0.5,
                                    star_pool_sat_off   = 0.5,
                                    star_pool_faint_off = 5.0,
                                    other_pix_frac_max  = 0.10,
                                    other_seg_pad       = 1L,
                                    bfgs_iters          = 200L,
                                    centre_tol_px       = 1.0,
                                    mag_tol_mag         = 0.1,
                                    chisq_max           = 3,
                                    do_2sig_clip        = TRUE,
                                    n_per_quadrant      = 2L,
                                    n_total_select      = 8L,
                                    galaxy_xy           = NULL,
                                    exclude_segids = NULL,
                                    exclude_galaxy_centre_seg = TRUE,
                                    target_flux_frac    = 0.99,
                                    default_fwhm_arcsec = NULL,
                                    default_con         = 2.1,
                                    subtract_segzero_sky = NULL,
                                    min_segzero_pix      = 50L,
                                    segim   = NULL,
                                    plot_profound = FALSE,
                                    ## ProFuse `profuseAllStarDoFit(..., plot=TRUE)` 4-panel montage
                                    ## (Data | Model | Data-Model | chi) on the v11-selected stars.
                                    plot_psf_fit = FALSE,
                                    optim_iters_psf_fit = 3L,
                                    Niters_psf_fit = c(200L, 200L),
                                    NfinalMCMC_psf_fit = 500L,
                                    seed_psf_fit = 666L,
                                    verbose = TRUE,
                                    output_kind = c("gama", "empirical",
                                                    "apodized", "default"),
                                    build_v9_sidecars = TRUE) {

  output_kind <- match.arg(output_kind)
  if (!is.null(segim)) stopifnot(identical(dim(image), dim(segim)))
  if (any(as.integer(psf_dim) %% 2L == 0L))
    psf_dim <- as.integer(psf_dim) + 1L
  if (!is.null(pixscale_arcsec)) {
    pixscale_arcsec <- as.numeric(pixscale_arcsec)[1]
    if (!is.finite(pixscale_arcsec) || pixscale_arcsec <= 0) pixscale_arcsec <- NULL
  }
  pixscale_for_attrs <- pixscale_arcsec %||% NA_real_

  ## ---- (a) global sky pedestal removal -----------------------------------
  sky_pedestal <- 0
  if (is.null(subtract_segzero_sky)) subtract_segzero_sky <- !is.null(segim)
  if (isTRUE(subtract_segzero_sky) && !is.null(segim)) {
    bg <- as.numeric(image[segim == 0L & is.finite(image)])
    if (length(bg) >= min_segzero_pix) {
      sky_pedestal <- stats::median(bg, na.rm = TRUE)
      if (is.finite(sky_pedestal) && abs(sky_pedestal) > 0) {
        image <- image - sky_pedestal
        if (verbose) message(sprintf(
          "[psf-v11] subtracted segim==0 sky pedestal = %.6g ADU", sky_pedestal))
      }
    }
  }

  ## ---- (b) ProFound segstats ---------------------------------------------
  pf_args <- list(image = image, magzero = magzero,
                  sky = 0, redosky = FALSE, boundstats = TRUE,
                  static_photom = TRUE, SBdilate = 2, skycut = 2,
                  verbose = FALSE, plot = plot_profound)
  if (!is.null(segim)) pf_args$segim <- segim
  pf <- do.call(ProFound::profoundProFound, pf_args)

  ## ---- (c) galaxy anchor & segIDs excluded from stellar pool --------------
  ## Galaxy / quadrant anchor (default = image centre).  Resolve *before* the
  ## star pool so we can drop `segim[centre]` — the science target is usually
  ## *not* a PSF star even if it slips the R50/R90 stellar locus.
  if (is.null(galaxy_xy))
    galaxy_xy <- c((nrow(image) + 1) / 2, (ncol(image) + 1) / 2)
  galaxy_centre_segID <- if (!is.null(segim))
    .segim_id_at_xy(segim, galaxy_xy) else NA_integer_
  star_pool_exclude_segids <- if (!is.null(exclude_segids)) {
    u <- unique(as.integer(exclude_segids[!is.na(exclude_segids)]))
    u[is.finite(u) & u > 0L]
  } else integer(0)
  if (isTRUE(exclude_galaxy_centre_seg) && is.finite(galaxy_centre_segID) &&
      galaxy_centre_segID > 0L)
    star_pool_exclude_segids <- unique(c(star_pool_exclude_segids, galaxy_centre_segID))
  if (verbose && length(star_pool_exclude_segids) > 0L) {
    msg <- sprintf("[psf-v11] PSF stellar pool excludes segID(s): %s",
                   paste(sort(star_pool_exclude_segids), collapse = ", "))
    if (isTRUE(exclude_galaxy_centre_seg) &&
        is.finite(galaxy_centre_segID) && galaxy_centre_segID > 0L &&
        galaxy_centre_segID %in% star_pool_exclude_segids)
      msg <- paste0(msg,
        sprintf(" (galaxy anchor xy=(%.2f,%.2f) -> segID %d)",
                galaxy_xy[1], galaxy_xy[2], galaxy_centre_segID))
    message(msg)
  }

  ## ---- (d) clean candidate pool with GAMA-style criteria -----------------
  cand <- select_clean_stars_v10(
    pf,
    n_max         = Nstar_pool,
    edge_excess_max = edge_excess_max,
    axrat_min     = star_pool_axrat_min,
    star_frac     = star_frac,
    sat_off       = star_pool_sat_off,
    faint_off     = star_pool_faint_off,
    exclude_segids = if (length(star_pool_exclude_segids))
      star_pool_exclude_segids else NULL,
    verbose       = verbose
  )

  ## ---- (e) per-star Moffat fits -----------------------------------------
  fits <- list()
  if (!is.null(cand) && nrow(cand) > 0L) {
    ## stamp half-size: ~6*FWHM from the segstats R50 (FWHM ~ 2*R50)
    R50_med <- stats::median(cand$R50, na.rm = TRUE)
    half_size <- max(8L, as.integer(round(6 * 2 * R50_med)))
    if (verbose) message(sprintf(
      "[psf-v11] per-star fit: half_size=%dpx (R50_med=%.2fpx)",
      half_size, R50_med))
    fits <- vector("list", nrow(cand))
    for (k in seq_len(nrow(cand))) {
      fits[[k]] <- .fit_one_star_gama_v11(
        image, segim, seg_row = cand[k, , drop = FALSE],
        half_size = half_size,
        other_seg_pad = other_seg_pad,
        bfgs_iters = bfgs_iters,
        magzero = magzero,
        other_pix_frac_max = other_pix_frac_max,
        verbose = verbose
      )
    }
  }

  ## ---- (f) suitability ---------------------------------------------------
  suitable <- if (length(fits) > 0L) {
    .gama_suitability(fits,
                      centre_tol_px = centre_tol_px,
                      mag_tol_mag   = mag_tol_mag,
                      chisq_max     = chisq_max,
                      do_2sig_clip  = do_2sig_clip,
                      verbose       = verbose)
  } else logical(0)

  ## ---- (g) per-star table ------------------------------------------------
  star_table <- if (length(fits) > 0L) {
    do.call(rbind, lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      pp <- if (isTRUE(f$ok)) f$fit_par else list()
      data.frame(
        idx = i,
        segID = f$segID %||% NA_integer_,
        ok = isTRUE(f$ok),
        suitable = if (length(suitable) >= i) suitable[i] else FALSE,
        reason = f$reason %||% NA_character_,
        segstats_xcen = f$segstats_xcen %||% NA_real_,
        segstats_ycen = f$segstats_ycen %||% NA_real_,
        segstats_mag  = f$segstats_mag  %||% NA_real_,
        segstats_R50  = f$segstats_R50  %||% NA_real_,
        segstats_R90  = f$segstats_R90  %||% NA_real_,
        fit_xcen = pp$xcen %||% NA_real_,
        fit_ycen = pp$ycen %||% NA_real_,
        fit_mag  = pp$mag  %||% NA_real_,
        fit_fwhm = pp$fwhm %||% NA_real_,
        fit_con  = pp$con  %||% NA_real_,
        fit_axrat = pp$axrat %||% NA_real_,
        fit_ang   = pp$ang  %||% NA_real_,
        fit_bg    = pp$bg   %||% NA_real_,
        reduced_chisq = f$reduced_chisq %||% NA_real_,
        ndof = f$ndof %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    }))
  } else NULL

  ## ---- (h) quadrant selection -------------------------------------------
  selected <- list(); selected_idx <- integer(0)
  if (length(fits) > 0L && any(suitable)) {
    suit_fits <- fits[suitable]
    suit_xy <- data.frame(
      xcen = vapply(suit_fits, function(f) f$segstats_xcen, numeric(1)),
      ycen = vapply(suit_fits, function(f) f$segstats_ycen, numeric(1))
    )
    pick <- .quadrant_select(suit_xy, galaxy_xy[1], galaxy_xy[2],
                             n_per_quad = n_per_quadrant,
                             n_total    = n_total_select)
    selected <- suit_fits[pick]
    ## map back to original index in `fits`
    suit_indices <- which(suitable)
    selected_idx <- suit_indices[pick]
  }

  ## ---- (i) build the GAMA median-parameter PSF ---------------------------
  psf_gama <- NULL
  median_seg_size_px <- NA_integer_
  if (length(selected) > 0L) {
    R50_sel <- vapply(selected, function(f) f$segstats_R50 %||% NA_real_,
                      numeric(1))
    if (any(is.finite(R50_sel)))
      median_seg_size_px <- as.integer(round(2 * stats::median(R50_sel,
                                                               na.rm = TRUE) * 6))
    psf_gama <- .build_median_param_psf(
      selected,
      target_flux_frac = target_flux_frac,
      max_dim = max(psf_dim),
      median_segment_size_px = median_seg_size_px,
      pixscale_arcsec = pixscale_for_attrs
    )
    if (verbose && !is.null(psf_gama))
      message(sprintf(
        "[psf-v11] GAMA model PSF: FWHM=%.3f px (%.3f arcsec)  beta=%.3f  axrat=%.3f  ang=%.1f  Nstar=%d  size=%dx%d",
        attr(psf_gama, "fwhm"), attr(psf_gama, "fwhm") *
          (pixscale_for_attrs %||% NA_real_),
        attr(psf_gama, "beta"), attr(psf_gama, "axrat"),
        attr(psf_gama, "ang"),  attr(psf_gama, "Nstar_used"),
        nrow(psf_gama), ncol(psf_gama)))
  }

  ## ---- (j) build the default round Moffat (always available as fallback) -
  default_fwhm_px <- if (!is.null(default_fwhm_arcsec) &&
                          !is.null(pixscale_arcsec)) {
    as.numeric(default_fwhm_arcsec) / pixscale_arcsec
  } else if (length(fits) > 0L) {
    fwhms_ok <- vapply(fits, function(f)
      if (isTRUE(f$ok)) f$fit_par$fwhm else NA_real_, numeric(1))
    if (any(is.finite(fwhms_ok))) stats::median(fwhms_ok, na.rm = TRUE)
    else 2.5
  } else 2.5
  psf_default <- .build_default_psf(default_fwhm_px, default_con,
                                    psf_dim = psf_dim,
                                    pixscale_arcsec = pixscale_for_attrs)

  ## ---- (k) PSF quality flag (GAMA bitwise) -------------------------------
  flag <- 0L
  if (is.null(cand) || nrow(cand) == 0L)               flag <- bitwOr(flag, 2L)
  if (length(fits) > 0L && !any(vapply(fits, function(f) isTRUE(f$ok),
                                       logical(1))))   flag <- bitwOr(flag, 4L)
  if (length(fits) > 0L && !any(suitable))             flag <- bitwOr(flag, 8L)
  if (length(selected) > 0L && length(selected) < n_total_select)
    flag <- bitwOr(flag, 16L)
  ## quadrant uneven: < n_per_quadrant in any quadrant
  if (length(selected) >= n_total_select) {
    sx <- vapply(selected, function(f) f$segstats_xcen, numeric(1))
    sy <- vapply(selected, function(f) f$segstats_ycen, numeric(1))
    qx <- ifelse(sx >= galaxy_xy[1], 1L, -1L)
    qy <- ifelse(sy >= galaxy_xy[2], 1L, -1L)
    qq <- table(paste(qx, qy))
    if (length(qq) < 4L || any(qq < n_per_quadrant))
      flag <- bitwOr(flag, 32L)
  }
  ## "default PSF used" bit
  if (is.null(psf_gama)) flag <- bitwOr(flag, 1L)
  ## DR4 PSFQUAL flag 64: median segment caps PSF size below 99% Moffat box
  if (!is.null(psf_gama) && !is.na(median_seg_size_px) &&
      median_seg_size_px > 0L) {
    fm <- attr(psf_gama, "fwhm")
    cn <- attr(psf_gama, "beta")
    if (is.finite(fm) && is.finite(cn)) {
      d_need <- .moffat_size_for_flux(fm, cn, target_flux_frac,
                                      max_dim = max(as.integer(psf_dim)),
                                      min_dim = 21L)
      if (is.finite(d_need) && d_need > median_seg_size_px)
        flag <- bitwOr(flag, 64L)
    }
  }

  ## ---- (l) sidecar v9 empirical / apodized for A/B comparison ------------
  psf_e <- NULL; stack_mof <- NULL; psf_apod <- NULL
  if (build_v9_sidecars && .has_v9() && length(selected) > 0L) {
    sel_xy <- data.frame(
      segID = vapply(selected, function(f) f$segID, integer(1)),
      xcen  = vapply(selected, function(f) f$segstats_xcen, numeric(1)),
      ycen  = vapply(selected, function(f) f$segstats_ycen, numeric(1)),
      mag   = vapply(selected, function(f) f$segstats_mag,  numeric(1)),
      R50   = vapply(selected, function(f) f$segstats_R50,  numeric(1)),
      axrat = vapply(selected, function(f) f$fit_par$axrat, numeric(1)),
      stringsAsFactors = FALSE
    )
    fwhm_for_stamp <- if (!is.null(psf_gama)) attr(psf_gama, "fwhm")
                      else 2.5
    psf_e <- tryCatch(
      build_empirical_psf_v9(image, segim, sel_xy,
                             psf_dim = psf_dim,
                             fwhm_px = fwhm_for_stamp,
                             refine_centroid = TRUE,
                             verbose = verbose),
      error = function(e) {
        if (verbose) message(sprintf(
          "[psf-v11] WARN sidecar empirical PSF failed: %s", conditionMessage(e)))
        NULL
      })
    if (!is.null(psf_e)) {
      stack_mof <- tryCatch(
        fit_moffat_to_stack(psf_e, allow_elliptical = TRUE,
                            fwhm_init = fwhm_for_stamp, beta_init = 3,
                            verbose = verbose),
        error = function(e) NULL)
      if (!is.null(stack_mof)) {
        attr(psf_e, "fwhm")  <- stack_mof$fwhm
        attr(psf_e, "beta")  <- stack_mof$beta
        attr(psf_e, "axrat") <- stack_mof$axrat
        attr(psf_e, "ang")   <- stack_mof$ang
        if (is.finite(pixscale_for_attrs))
          attr(psf_e, "fwhm_arcsec") <- stack_mof$fwhm * pixscale_for_attrs
        psf_apod <- tryCatch(
          build_apodized_psf(psf_e, stack_mof,
                             r1_fwhm = 3.0, r2_fwhm = 5.0, verbose = verbose),
          error = function(e) NULL)
        if (!is.null(psf_apod) && is.finite(pixscale_for_attrs)) {
          attr(psf_apod, "pixscale_arcsec") <- pixscale_for_attrs
          attr(psf_apod, "fwhm_arcsec") <- attr(psf_apod, "fwhm") * pixscale_for_attrs
        }
      }
    }
  }

  ## ---- (m) recommended PSF -----------------------------------------------
  recommended <- switch(
    output_kind,
    "gama"      = psf_gama %||% psf_default,
    "empirical" = psf_e %||% psf_gama %||% psf_default,
    "apodized"  = psf_apod %||% psf_gama %||% psf_default,
    "default"   = psf_default
  )
  output_kind_used <- if (identical(recommended, psf_gama)    && !is.null(psf_gama))    "gama"
                      else if (identical(recommended, psf_e)  && !is.null(psf_e))      "empirical"
                      else if (identical(recommended, psf_apod) && !is.null(psf_apod)) "apodized"
                      else                                                              "default"
  if (verbose)
    message(sprintf(
      "[psf-v11] output=%s | gama=%s | empirical=%s | default FWHM=%.3fpx (con=%.2f) | flag=%d",
      output_kind_used,
      if (is.null(psf_gama)) "n/a" else sprintf("FWHM=%.3f beta=%.3f Nstar=%d",
        attr(psf_gama, "fwhm"), attr(psf_gama, "beta"), attr(psf_gama, "Nstar_used")),
      if (is.null(psf_e)) "n/a" else sprintf("FWHM=%.3f Nstar=%d",
        attr(psf_e, "fwhm"), attr(psf_e, "Nstar_used") %||% NA_integer_),
      default_fwhm_px, default_con, flag))

  ## ProFuse all-star notebook diagnostic (build_psfs v9: PLOT_PSF_FIT ->
  ## plot_psf_fit). Quadrant-selected stars first; if none yet the pool /
  ## converged fits still allows the default ProFuse Data|Model|resid montage.
  if (isTRUE(plot_psf_fit)) {
    locs <- NULL
    star_dom_mag <- NA_real_
    star_circ <- TRUE

    if (!is.null(star_table) && length(selected_idx) > 0L) {
      sel_vis <- star_table[star_table$idx %in% selected_idx, , drop = FALSE]
      if (nrow(sel_vis) > 0L) {
        star_dom_mag <- max(sel_vis$segstats_mag, na.rm = TRUE) + 0.01
        locs <- cbind(sel_vis$segstats_xcen, sel_vis$segstats_ycen)
        ax_med <- stats::median(sel_vis$fit_axrat, na.rm = TRUE)
        star_circ <- !is.finite(ax_med) || ax_med >= 0.95
      }
    }

    if (is.null(locs) || nrow(locs) == 0L) {
      if (!is.null(star_table) && nrow(star_table) > 0L &&
          "ok" %in% names(star_table)) {
        ok_tab <- star_table[which(star_table$ok %in% TRUE), , drop = FALSE]
        if (nrow(ok_tab) > 0L) {
          ok_tab <- ok_tab[order(ok_tab$segstats_mag), , drop = FALSE]
          nloc <- min(as.integer(n_total_select)[1], nrow(ok_tab))
          ok_tab <- ok_tab[seq_len(nloc), , drop = FALSE]
          star_dom_mag <- max(ok_tab$segstats_mag, na.rm = TRUE) + 0.01
          locs <- cbind(ok_tab$segstats_xcen, ok_tab$segstats_ycen)
          ax_med <- stats::median(ok_tab$fit_axrat, na.rm = TRUE)
          star_circ <- !is.finite(ax_med) || ax_med >= 0.95
        }
      }
    }

    if (is.null(locs) || nrow(locs) == 0L) {
      if (!is.null(cand) && nrow(cand) > 0L &&
          all(c("xcen", "ycen", "mag") %in% names(cand))) {
        cand_s <- cand[order(cand$mag), , drop = FALSE]
        nloc <- min(as.integer(n_total_select)[1], nrow(cand_s))
        cand_s <- cand_s[seq_len(nloc), , drop = FALSE]
        star_dom_mag <- max(cand_s$mag, na.rm = TRUE) + 0.01
        locs <- cbind(as.numeric(cand_s$xcen), as.numeric(cand_s$ycen))
        if ("axrat" %in% names(cand_s)) {
          ax_med <- stats::median(as.numeric(cand_s$axrat), na.rm = TRUE)
          star_circ <- !is.finite(ax_med) || ax_med >= 0.95
        } else star_circ <- TRUE
      }
    }

    if (!is.null(locs) && nrow(locs) > 0L && is.finite(star_dom_mag) &&
        all(is.finite(locs))) {
      colnames(locs) <- NULL
      allstar_viz_args <- list(
        image         = image,
        magzero       = magzero,
        psf_dim       = psf_dim,
        Nstar         = nrow(locs),
        locs          = locs,
        star_dom_mag  = star_dom_mag,
        star_con      = 2.5,
        star_con_fit  = TRUE,
        star_circ     = star_circ,
        SBdilate      = 2,
        rough         = FALSE,
        autoclip      = is.null(segim),
        roughpedestal = TRUE,
        optim_iters   = as.integer(optim_iters_psf_fit)[1],
        Niters        = as.integer(Niters_psf_fit),
        NfinalMCMC    = as.integer(NfinalMCMC_psf_fit)[1],
        seed          = as.integer(seed_psf_fit)[1],
        plot          = TRUE
      )
      if (!is.null(segim)) allstar_viz_args$segim <- segim
      tryCatch(
        {
          do.call(ProFuse::profuseAllStarDoFit, allstar_viz_args)
          if (verbose)
            message("[psf-v11] ProFuse all-star montage plotted (profuseAllStarDoFit plot=TRUE)")
        },
        error = function(e) {
          if (verbose)
            message(sprintf(
              "[psf-v11] profuseAllStarDoFit (plot) failed: %s",
              conditionMessage(e)))
        }
      )
    }
  }

  list(
    psf                    = recommended,
    psf_gama               = psf_gama,
    psf_default            = psf_default,
    psf_per_star_models    = lapply(selected, function(f) {
      pp <- f$fit_par
      list(fwhm = pp$fwhm, con = pp$con, axrat = pp$axrat, ang = pp$ang,
           xcen = pp$xcen, ycen = pp$ycen, mag = pp$mag,
           reduced_chisq = f$reduced_chisq, segID = f$segID)
    }),
    psf_empirical          = psf_e,
    psf_apodized           = psf_apod,
    stack_moffat           = stack_mof,
    output_kind            = output_kind_used,
    profound               = pf,
    candidates             = cand,
    fits                   = fits,
    suitable               = suitable,
    selected_idx           = selected_idx,
    star_table             = star_table,
    decision               = list(
      Nstar_pool         = Nstar_pool,
      star_frac          = star_frac,
      edge_excess_max    = edge_excess_max,
      axrat_min          = star_pool_axrat_min,
      sat_off            = star_pool_sat_off,
      faint_off          = star_pool_faint_off,
      other_pix_frac_max = other_pix_frac_max,
      bfgs_iters         = bfgs_iters,
      centre_tol_px      = centre_tol_px,
      mag_tol_mag        = mag_tol_mag,
      chisq_max          = chisq_max,
      do_2sig_clip       = do_2sig_clip,
      n_per_quadrant     = n_per_quadrant,
      n_total_select     = n_total_select,
      target_flux_frac   = target_flux_frac,
      default_fwhm_px    = default_fwhm_px,
      default_con        = default_con,
      median_seg_size_px = median_seg_size_px,
      galaxy_xy          = galaxy_xy,
      star_pool_exclude_segids = star_pool_exclude_segids,
      galaxy_centre_segID      = galaxy_centre_segID,
      exclude_galaxy_centre_seg = exclude_galaxy_centre_seg
    ),
    psf_quality_flag        = as.integer(flag),
    segim                   = segim,
    stars                   = if (!is.null(star_table)) {
      sel <- star_table[star_table$idx %in% selected_idx, , drop = FALSE]
      data.frame(
        segID = sel$segID,
        xcen  = sel$segstats_xcen,
        ycen  = sel$segstats_ycen,
        mag   = sel$segstats_mag,
        R50   = sel$segstats_R50,
        axrat = sel$fit_axrat,
        stringsAsFactors = FALSE
      )
    } else data.frame(),
    sky_pedestal_subtracted = as.numeric(sky_pedestal)
  )
}

## ----------------------------------------------------------------------------
##  Inline 2×4 PSF montage for IRkernel (drawn *before* profuseAllStarDoFit so it
##  is not buried after dozens of MCMC / chi panels).
## ----------------------------------------------------------------------------
.notebook_psf_magimage_2x4 <- function(psf_gama, psf_default, psf_e, psf_apod,
                                     pixscale_arcsec = NA_real_) {
  pl <- function(p, lab) {
    if (is.null(p)) {
      graphics::plot.new()
      graphics::title(main = sprintf("%s (n/a)", lab))
      return(invisible())
    }
    fwhm <- attr(p, "fwhm") %||% NA_real_
    beta <- attr(p, "beta") %||% NA_real_
    axr  <- attr(p, "axrat") %||% NA_real_
    sub  <- sprintf("FWHM=%.2fpx  beta=%.2f  axrat=%.2f", fwhm, beta, axr)
    if (is.finite(pixscale_arcsec) && is.finite(fwhm))
      sub <- paste0(sub, sprintf("  (%.2f arcsec)", fwhm * pixscale_arcsec))
    magicaxis::magimage(p, main = sprintf("%s linear", lab), sub = sub)
    shift <- max(1e-12, -min(p, na.rm = TRUE) + 1e-12)
    magicaxis::magimage(log10(p + shift), main = sprintf("%s log", lab), sub = sub)
  }
  op <- graphics::par(mfrow = c(2, 4), mar = c(2, 2, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  pl(psf_gama,    "GAMA median-param Moffat")
  pl(psf_default, "default round Moffat")
  pl(psf_e,       "v9 empirical (recentred)")
  pl(psf_apod,    "v9 apodized hybrid")
  invisible(TRUE)
}

## ============================================================================
##  build_psf_for_image_v13 — BDInputsRv05-style path:
##    * strict stellar pool (`select_clean_stars_v10` once; gri four percent + skycut=2);
##    * `.gama_suitability_dr4` — 2sigma on FWHM, con−1, axrat, angle (no bg);
##    * `profit_optim_dr4` selects BFGS (GAMA) vs L-BFGS-B (WAVES default).
## ============================================================================
build_psf_for_image_v13 <- function(image, magzero,
                                    pixscale_arcsec = NULL,
                                    psf_dim = c(75L, 75L),
                                    Nstar_pool = 50L,
                                    star_frac = 0.04,
                                    edge_excess_max = 1.0,
                                    star_pool_axrat_min = 0.5,
                                    star_pool_sat_off   = 0.5,
                                    star_pool_faint_off = 5.0,
                                    other_pix_frac_max  = 0.10,
                                    other_seg_pad       = 1L,
                                    bfgs_iters          = 200L,
                                    centre_tol_px       = 1.0,
                                    mag_tol_mag         = 0.1,
                                    chisq_max           = 3,
                                    do_2sig_clip        = TRUE,
                                    n_per_quadrant      = 2L,
                                    n_total_select      = 8L,
                                    galaxy_xy           = NULL,
                                    exclude_segids = NULL,
                                    exclude_galaxy_centre_seg = TRUE,
                                    target_flux_frac    = 0.99,
                                    default_fwhm_arcsec = NULL,
                                    default_con         = 2.1,
                                    subtract_segzero_sky = NULL,
                                    min_segzero_pix      = 50L,
                                    segim   = NULL,
                                    plot_profound = FALSE,
                                    plot_psf_fit = FALSE,
                                    optim_iters_psf_fit = 3L,
                                    Niters_psf_fit = c(200L, 200L),
                                    NfinalMCMC_psf_fit = 500L,
                                    seed_psf_fit = 666L,
                                    verbose = TRUE,
                                    output_kind = c("gama", "empirical",
                                                    "apodized", "default"),
                                    profit_optim_dr4 = c("L-BFGS-B", "BFGS"),
                                    psfqual_skycut_segstats_bounds = NULL,
                                    build_v9_sidecars = TRUE) {

  output_kind <- match.arg(output_kind)
  profit_optim_dr4 <- match.arg(profit_optim_dr4)
  fit_star_fn <- if (identical(profit_optim_dr4, "BFGS"))
    .fit_one_star_gama_v11 else .fit_one_star_gama_v12
  if (!is.null(segim)) stopifnot(identical(dim(image), dim(segim)))
  if (any(as.integer(psf_dim) %% 2L == 0L))
    psf_dim <- as.integer(psf_dim) + 1L
  if (!is.null(pixscale_arcsec)) {
    pixscale_arcsec <- as.numeric(pixscale_arcsec)[1]
    if (!is.finite(pixscale_arcsec) || pixscale_arcsec <= 0) pixscale_arcsec <- NULL
  }
  pixscale_for_attrs <- pixscale_arcsec %||% NA_real_

  ## (a) global sky pedestal removal (same as v11)
  sky_pedestal <- 0
  if (is.null(subtract_segzero_sky)) subtract_segzero_sky <- !is.null(segim)
  if (isTRUE(subtract_segzero_sky) && !is.null(segim)) {
    bg <- as.numeric(image[segim == 0L & is.finite(image)])
    if (length(bg) >= min_segzero_pix) {
      sky_pedestal <- stats::median(bg, na.rm = TRUE)
      if (is.finite(sky_pedestal) && abs(sky_pedestal) > 0) {
        image <- image - sky_pedestal
        if (verbose) message(sprintf(
          "[psf-v13] subtracted segim==0 sky pedestal = %.6g ADU", sky_pedestal))
      }
    }
  }

  ## (b) ProFound segstats
  pf_args <- list(image = image, magzero = magzero,
                  sky = 0, redosky = FALSE, boundstats = TRUE,
                  static_photom = TRUE, SBdilate = 2, skycut = 2,
                  verbose = FALSE, plot = plot_profound)
  if (!is.null(segim)) pf_args$segim <- segim
  pf <- do.call(ProFound::profoundProFound, pf_args)

  n_segstats_pf <- if (!is.null(pf$segstats) && nrow(pf$segstats) > 0L)
    nrow(pf$segstats) else 0L

  ## (c) Galaxy anchor & segIDs excluded from stellar pool
  if (is.null(galaxy_xy))
    galaxy_xy <- c((nrow(image) + 1) / 2, (ncol(image) + 1) / 2)
  galaxy_centre_segID <- if (!is.null(segim))
    .segim_id_at_xy(segim, galaxy_xy) else NA_integer_
  star_pool_exclude_segids <- if (!is.null(exclude_segids)) {
    u <- unique(as.integer(exclude_segids[!is.na(exclude_segids)]))
    u[is.finite(u) & u > 0L]
  } else integer(0)
  if (isTRUE(exclude_galaxy_centre_seg) && is.finite(galaxy_centre_segID) &&
      galaxy_centre_segID > 0L)
    star_pool_exclude_segids <- unique(c(star_pool_exclude_segids, galaxy_centre_segID))
  if (verbose && length(star_pool_exclude_segids) > 0L) {
    msg <- sprintf("[psf-v13] PSF stellar pool excludes segID(s): %s",
                   paste(sort(star_pool_exclude_segids), collapse = ", "))
    if (isTRUE(exclude_galaxy_centre_seg) &&
        is.finite(galaxy_centre_segID) && galaxy_centre_segID > 0L &&
        galaxy_centre_segID %in% star_pool_exclude_segids)
      msg <- paste0(msg,
        sprintf(" (galaxy anchor xy=(%.2f,%.2f) -> segID %d)",
                galaxy_xy[1], galaxy_xy[2], galaxy_centre_segID))
    message(msg)
  }

  ## (d) Strict BDInputsRv05 stellar pool (select_clean_stars_v10 once; gri ~4%)
  cand <- select_clean_stars_v10(
    pf,
    n_max         = Nstar_pool,
    edge_excess_max = edge_excess_max,
    axrat_min     = star_pool_axrat_min,
    star_frac     = star_frac,
    sat_off       = star_pool_sat_off,
    faint_off     = star_pool_faint_off,
    exclude_segids = if (length(star_pool_exclude_segids))
      star_pool_exclude_segids else NULL,
    verbose       = verbose
  )
  pool_tier <- NA_integer_
  pool_tier_label <- "BDInputsRv05_strict_v10_only"

  ## (e) per-star Moffat (ProFit; GAMA=BFGS, WAVES default L-BFGS-B)
  fits <- list()
  if (!is.null(cand) && nrow(cand) > 0L) {
    R50_med <- stats::median(cand$R50, na.rm = TRUE)
    half_size <- max(8L, as.integer(round(6 * 2 * R50_med)))
    if (verbose) message(sprintf(
      "[psf-v13] per-star fit: half_size=%dpx (R50_med=%.2fpx)",
      half_size, R50_med))
    fits <- vector("list", nrow(cand))
    for (k in seq_len(nrow(cand))) {
      fits[[k]] <- fit_star_fn(
        image, segim, seg_row = cand[k, , drop = FALSE],
        half_size = half_size,
        other_seg_pad = other_seg_pad,
        bfgs_iters = bfgs_iters,
        magzero = magzero,
        other_pix_frac_max = other_pix_frac_max,
        verbose = verbose
      )
    }
  }

  ## (f) suitability (BDInputs post-fit gates + 2-sigma clip w/o bg)
  suitable <- if (length(fits) > 0L) {
    .gama_suitability_dr4(fits,
                      centre_tol_px = centre_tol_px,
                      mag_tol_mag   = mag_tol_mag,
                      chisq_max     = chisq_max,
                      do_2sig_clip  = do_2sig_clip,
                      verbose       = verbose)
  } else logical(0)

  ## (g) per-star table (same as v11)
  star_table <- if (length(fits) > 0L) {
    do.call(rbind, lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      pp <- if (isTRUE(f$ok)) f$fit_par else list()
      data.frame(
        idx = i,
        segID = f$segID %||% NA_integer_,
        ok = isTRUE(f$ok),
        suitable = if (length(suitable) >= i) suitable[i] else FALSE,
        reason = f$reason %||% NA_character_,
        segstats_xcen = f$segstats_xcen %||% NA_real_,
        segstats_ycen = f$segstats_ycen %||% NA_real_,
        segstats_mag  = f$segstats_mag  %||% NA_real_,
        segstats_R50  = f$segstats_R50  %||% NA_real_,
        segstats_R90  = f$segstats_R90  %||% NA_real_,
        fit_xcen = pp$xcen %||% NA_real_,
        fit_ycen = pp$ycen %||% NA_real_,
        fit_mag  = pp$mag  %||% NA_real_,
        fit_fwhm = pp$fwhm %||% NA_real_,
        fit_con  = pp$con  %||% NA_real_,
        fit_axrat = pp$axrat %||% NA_real_,
        fit_ang   = pp$ang  %||% NA_real_,
        fit_bg    = pp$bg   %||% NA_real_,
        reduced_chisq = f$reduced_chisq %||% NA_real_,
        ndof = f$ndof %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    }))
  } else NULL

  ## (h) quadrant selection (same as v11)
  selected <- list(); selected_idx <- integer(0)
  if (length(fits) > 0L && any(suitable)) {
    suit_fits <- fits[suitable]
    suit_xy <- data.frame(
      xcen = vapply(suit_fits, function(f) f$segstats_xcen, numeric(1)),
      ycen = vapply(suit_fits, function(f) f$segstats_ycen, numeric(1))
    )
    pick <- .quadrant_select(suit_xy, galaxy_xy[1], galaxy_xy[2],
                             n_per_quad = n_per_quadrant,
                             n_total    = n_total_select)
    selected <- suit_fits[pick]
    suit_indices <- which(suitable)
    selected_idx <- suit_indices[pick]
  }

  ## (i) build the GAMA median-parameter PSF (same as v11)
  psf_gama <- NULL
  median_seg_size_px <- NA_integer_
  if (length(selected) > 0L) {
    R50_sel <- vapply(selected, function(f) f$segstats_R50 %||% NA_real_,
                      numeric(1))
    if (any(is.finite(R50_sel)))
      median_seg_size_px <- as.integer(round(2 * stats::median(R50_sel,
                                                               na.rm = TRUE) * 6))
    psf_gama <- .build_median_param_psf(
      selected,
      target_flux_frac = target_flux_frac,
      max_dim = max(psf_dim),
      median_segment_size_px = median_seg_size_px,
      pixscale_arcsec = pixscale_for_attrs
    )
    if (verbose && !is.null(psf_gama))
      message(sprintf(
        "[psf-v13] GAMA model PSF: FWHM=%.3f px (%.3f arcsec)  beta=%.3f  axrat=%.3f  ang=%.1f  Nstar=%d  size=%dx%d",
        attr(psf_gama, "fwhm"),
        attr(psf_gama, "fwhm") * (pixscale_for_attrs %||% NA_real_),
        attr(psf_gama, "beta"), attr(psf_gama, "axrat"),
        attr(psf_gama, "ang"),  attr(psf_gama, "Nstar_used"),
        nrow(psf_gama), ncol(psf_gama)))
  }

  ## (j) default round Moffat
  default_fwhm_px <- if (!is.null(default_fwhm_arcsec) &&
                          !is.null(pixscale_arcsec)) {
    as.numeric(default_fwhm_arcsec) / pixscale_arcsec
  } else if (length(fits) > 0L) {
    fwhms_ok <- vapply(fits, function(f)
      if (isTRUE(f$ok)) f$fit_par$fwhm else NA_real_, numeric(1))
    if (any(is.finite(fwhms_ok))) stats::median(fwhms_ok, na.rm = TRUE)
    else 2.5
  } else 2.5
  psf_default <- .build_default_psf(default_fwhm_px, default_con,
                                    psf_dim = psf_dim,
                                    pixscale_arcsec = pixscale_for_attrs)

  ## (k) DR4 PSFQUAL_FLAG (BDInputsRv05)
  flag <- 0L
  ## 2 = no sources in ProFound segstats; 4 = sources but none survive star pool
  if (n_segstats_pf == 0L)                               flag <- bitwOr(flag, 2L)
  else if (is.null(cand) || nrow(cand) == 0L)           flag <- bitwOr(flag, 4L)
  if (length(fits) > 0L && !any(vapply(fits, function(f) isTRUE(f$ok),
                                       logical(1))))   flag <- bitwOr(flag, 4L)
  if (length(fits) > 0L && !any(suitable))             flag <- bitwOr(flag, 8L)
  if (length(selected) > 0L && length(selected) < n_total_select)
    flag <- bitwOr(flag, 16L)
  if (length(selected) >= n_total_select) {
    sx <- vapply(selected, function(f) f$segstats_xcen, numeric(1))
    sy <- vapply(selected, function(f) f$segstats_ycen, numeric(1))
    qx <- ifelse(sx >= galaxy_xy[1], 1L, -1L)
    qy <- ifelse(sy >= galaxy_xy[2], 1L, -1L)
    qq <- table(paste(qx, qy))
    if (length(qq) < 4L || any(qq < n_per_quadrant))
      flag <- bitwOr(flag, 32L)
  }
  if (is.null(psf_gama)) flag <- bitwOr(flag, 1L)
  if (!is.null(psf_gama) && !is.na(median_seg_size_px) &&
      median_seg_size_px > 0L) {
    fm <- attr(psf_gama, "fwhm")
    cn <- attr(psf_gama, "beta")
    if (is.finite(fm) && is.finite(cn)) {
      d_need <- .moffat_size_for_flux(fm, cn, target_flux_frac,
                                      max_dim = max(as.integer(psf_dim)),
                                      min_dim = 21L)
      if (is.finite(d_need) && d_need > median_seg_size_px)
        flag <- bitwOr(flag, 64L)
    }
  }
  ## DR4 PSFQUAL bit 128: optional skycut sanity from ProFound segstats count
  if (!is.null(psfqual_skycut_segstats_bounds) &&
      length(psfqual_skycut_segstats_bounds) >= 2L) {
    lo <- suppressWarnings(as.integer(psfqual_skycut_segstats_bounds[1]))
    hi <- suppressWarnings(as.integer(psfqual_skycut_segstats_bounds[2]))
    if (is.finite(lo) && is.finite(hi) && hi >= lo &&
        (n_segstats_pf < lo || n_segstats_pf > hi))
      flag <- bitwOr(flag, 128L)
  }

  ## (l) sidecar v9 empirical / apodized (same as v11)
  psf_e <- NULL; stack_mof <- NULL; psf_apod <- NULL
  if (build_v9_sidecars && .has_v9() && length(selected) > 0L) {
    sel_xy <- data.frame(
      segID = vapply(selected, function(f) f$segID, integer(1)),
      xcen  = vapply(selected, function(f) f$segstats_xcen, numeric(1)),
      ycen  = vapply(selected, function(f) f$segstats_ycen, numeric(1)),
      mag   = vapply(selected, function(f) f$segstats_mag,  numeric(1)),
      R50   = vapply(selected, function(f) f$segstats_R50,  numeric(1)),
      axrat = vapply(selected, function(f) f$fit_par$axrat, numeric(1)),
      stringsAsFactors = FALSE
    )
    fwhm_for_stamp <- if (!is.null(psf_gama)) attr(psf_gama, "fwhm")
                      else 2.5
    psf_e <- tryCatch(
      build_empirical_psf_v9(image, segim, sel_xy,
                             psf_dim = psf_dim,
                             fwhm_px = fwhm_for_stamp,
                             refine_centroid = TRUE,
                             verbose = verbose),
      error = function(e) {
        if (verbose) message(sprintf(
          "[psf-v13] WARN sidecar empirical PSF failed: %s", conditionMessage(e)))
        NULL
      })
    if (!is.null(psf_e)) {
      stack_mof <- tryCatch(
        fit_moffat_to_stack(psf_e, allow_elliptical = TRUE,
                            fwhm_init = fwhm_for_stamp, beta_init = 3,
                            verbose = verbose),
        error = function(e) NULL)
      if (!is.null(stack_mof)) {
        attr(psf_e, "fwhm")  <- stack_mof$fwhm
        attr(psf_e, "beta")  <- stack_mof$beta
        attr(psf_e, "axrat") <- stack_mof$axrat
        attr(psf_e, "ang")   <- stack_mof$ang
        if (is.finite(pixscale_for_attrs))
          attr(psf_e, "fwhm_arcsec") <- stack_mof$fwhm * pixscale_for_attrs
        psf_apod <- tryCatch(
          build_apodized_psf(psf_e, stack_mof,
                             r1_fwhm = 3.0, r2_fwhm = 5.0, verbose = verbose),
          error = function(e) NULL)
        if (!is.null(psf_apod) && is.finite(pixscale_for_attrs)) {
          attr(psf_apod, "pixscale_arcsec") <- pixscale_for_attrs
          attr(psf_apod, "fwhm_arcsec") <- attr(psf_apod, "fwhm") * pixscale_for_attrs
        }
      }
    }
  }

  ## (m) recommended PSF (same as v11)
  recommended <- switch(
    output_kind,
    "gama"      = psf_gama %||% psf_default,
    "empirical" = psf_e %||% psf_gama %||% psf_default,
    "apodized"  = psf_apod %||% psf_gama %||% psf_default,
    "default"   = psf_default
  )
  output_kind_used <- if (identical(recommended, psf_gama)    && !is.null(psf_gama))    "gama"
                      else if (identical(recommended, psf_e)  && !is.null(psf_e))      "empirical"
                      else if (identical(recommended, psf_apod) && !is.null(psf_apod)) "apodized"
                      else                                                              "default"
  if (verbose)
    message(sprintf(
      "[psf-v13] output=%s | gama=%s | empirical=%s | default FWHM=%.3fpx (con=%.2f) | flag=%d | optim=%s | Nsegstats=%d",
      output_kind_used,
      if (is.null(psf_gama)) "n/a" else sprintf("FWHM=%.3f beta=%.3f Nstar=%d",
        attr(psf_gama, "fwhm"), attr(psf_gama, "beta"), attr(psf_gama, "Nstar_used")),
      if (is.null(psf_e)) "n/a" else sprintf("FWHM=%.3f Nstar=%d",
        attr(psf_e, "fwhm"), attr(psf_e, "Nstar_used") %||% NA_integer_),
      default_fwhm_px, default_con, flag, profit_optim_dr4, n_segstats_pf))

  ## PSF montage (notebook figure) before heavy ProFuse data/model/MCMC plots
  psf_montage_drawn_early <- FALSE
  if (isTRUE(plot_psf_fit) && requireNamespace("magicaxis", quietly = TRUE)) {
    tryCatch({
      .notebook_psf_magimage_2x4(
        psf_gama, psf_default, psf_e, psf_apod,
        pixscale_arcsec = pixscale_for_attrs %||% NA_real_)
      psf_montage_drawn_early <- TRUE
      if (verbose)
        message("[psf-v13] notebook PSF montage drawn before profuseAllStarDoFit.")
    }, error = function(e) invisible(NULL))
  }

  ## ProFuse all-star diagnostic (Casura PSF appendix style). Same fallbacks as v11.
  if (isTRUE(plot_psf_fit)) {
    locs <- NULL
    star_dom_mag <- NA_real_
    star_circ <- TRUE

    if (!is.null(star_table) && length(selected_idx) > 0L) {
      sel_vis <- star_table[star_table$idx %in% selected_idx, , drop = FALSE]
      if (nrow(sel_vis) > 0L) {
        star_dom_mag <- max(sel_vis$segstats_mag, na.rm = TRUE) + 0.01
        locs <- cbind(sel_vis$segstats_xcen, sel_vis$segstats_ycen)
        ax_med <- stats::median(sel_vis$fit_axrat, na.rm = TRUE)
        star_circ <- !is.finite(ax_med) || ax_med >= 0.95
      }
    }

    if (is.null(locs) || nrow(locs) == 0L) {
      if (!is.null(star_table) && nrow(star_table) > 0L &&
          "ok" %in% names(star_table)) {
        ok_tab <- star_table[which(star_table$ok %in% TRUE), , drop = FALSE]
        if (nrow(ok_tab) > 0L) {
          ok_tab <- ok_tab[order(ok_tab$segstats_mag), , drop = FALSE]
          nloc <- min(as.integer(n_total_select)[1], nrow(ok_tab))
          ok_tab <- ok_tab[seq_len(nloc), , drop = FALSE]
          star_dom_mag <- max(ok_tab$segstats_mag, na.rm = TRUE) + 0.01
          locs <- cbind(ok_tab$segstats_xcen, ok_tab$segstats_ycen)
          ax_med <- stats::median(ok_tab$fit_axrat, na.rm = TRUE)
          star_circ <- !is.finite(ax_med) || ax_med >= 0.95
        }
      }
    }

    if (is.null(locs) || nrow(locs) == 0L) {
      if (!is.null(cand) && nrow(cand) > 0L &&
          all(c("xcen", "ycen", "mag") %in% names(cand))) {
        cand_s <- cand[order(cand$mag), , drop = FALSE]
        nloc <- min(as.integer(n_total_select)[1], nrow(cand_s))
        cand_s <- cand_s[seq_len(nloc), , drop = FALSE]
        star_dom_mag <- max(cand_s$mag, na.rm = TRUE) + 0.01
        locs <- cbind(as.numeric(cand_s$xcen), as.numeric(cand_s$ycen))
        if ("axrat" %in% names(cand_s)) {
          ax_med <- stats::median(as.numeric(cand_s$axrat), na.rm = TRUE)
          star_circ <- !is.finite(ax_med) || ax_med >= 0.95
        } else {
          star_circ <- TRUE
        }
      }
    }

    pf_exists <- suppressWarnings(requireNamespace("ProFuse", quietly = TRUE)) &&
                  exists("profuseAllStarDoFit",
                         envir = asNamespace("ProFuse"), inherits = FALSE)

    if (!pf_exists && verbose) {
      message("[psf-v13] ProFuse not available; skipping profuseAllStarDoFit.")
    }

    if (pf_exists && !is.null(locs) && nrow(locs) > 0L && is.finite(star_dom_mag) &&
        all(is.finite(locs))) {
      colnames(locs) <- NULL
      Niters_pf <- Niters_psf_fit
      if (!is.finite(Niters_pf[1]) || length(Niters_pf) < 2L) {
        Niters_pf <- rep(200L, 2L)
      } else Niters_pf <- as.integer(Niters_pf)[1:2]

      allstar_viz_args <- list(
        image         = image,
        magzero       = magzero,
        psf_dim       = psf_dim,
        Nstar         = nrow(locs),
        locs          = locs,
        star_dom_mag  = star_dom_mag,
        star_con      = 2.5,
        star_con_fit  = TRUE,
        star_circ     = star_circ,
        SBdilate      = 2,
        rough         = FALSE,
        autoclip      = is.null(segim),
        roughpedestal = TRUE,
        optim_iters   = as.integer(optim_iters_psf_fit)[1],
        Niters        = Niters_pf,
        NfinalMCMC    = as.integer(NfinalMCMC_psf_fit)[1],
        seed          = as.integer(seed_psf_fit)[1],
        plot          = TRUE)
      if (!is.null(segim)) allstar_viz_args$segim <- segim
      tryCatch({
          do.call(ProFuse::profuseAllStarDoFit, allstar_viz_args)
          if (verbose)
            message("[psf-v13] profuseAllStarDoFit plotted (plot=TRUE)")
        }, error = function(e) {
          if (verbose)
            message(sprintf("[psf-v13] profuseAllStarDoFit failed: %s",
                           conditionMessage(e)))
        })
    }
  }

  list(
    psf                    = recommended,
    psf_gama               = psf_gama,
    psf_default            = psf_default,
    psf_per_star_models    = lapply(selected, function(f) {
      pp <- f$fit_par
      list(
        fwhm = pp$fwhm, con = pp$con, axrat = pp$axrat, ang = pp$ang,
        xcen = pp$xcen, ycen = pp$ycen, mag = pp$mag,
        reduced_chisq = f$reduced_chisq, segID = f$segID)
    }),
    psf_empirical          = psf_e,
    psf_apodized           = psf_apod,
    stack_moffat           = stack_mof,
    output_kind            = output_kind_used,
    profound               = pf,
    candidates             = cand,
    fits                   = fits,
    suitable               = suitable,
    selected_idx           = selected_idx,
    star_table             = star_table,
    decision               = list(
      Nstar_pool         = Nstar_pool,
      pool_tier          = pool_tier,
      pool_tier_label    = pool_tier_label,
      star_frac          = star_frac,
      edge_excess_max    = edge_excess_max,
      axrat_min          = star_pool_axrat_min,
      sat_off            = star_pool_sat_off,
      faint_off          = star_pool_faint_off,
      other_pix_frac_max = other_pix_frac_max,
      bfgs_iters         = bfgs_iters,
      centre_tol_px      = centre_tol_px,
      mag_tol_mag        = mag_tol_mag,
      chisq_max          = chisq_max,
      do_2sig_clip       = do_2sig_clip,
      n_per_quadrant     = n_per_quadrant,
      n_total_select     = n_total_select,
      target_flux_frac   = target_flux_frac,
      default_fwhm_px    = default_fwhm_px,
      default_con        = default_con,
      median_seg_size_px = median_seg_size_px,
      galaxy_xy          = galaxy_xy,
      star_pool_exclude_segids = star_pool_exclude_segids,
      galaxy_centre_segID      = galaxy_centre_segID,
      exclude_galaxy_centre_seg = exclude_galaxy_centre_seg,
      profit_optim_dr4          = profit_optim_dr4,
      profound_segstats_n       = n_segstats_pf,
      psfqual_skycut_segstats_bounds = psfqual_skycut_segstats_bounds),
    psf_quality_flag       = as.integer(flag),
    segim                   = segim,
    stars                   = if (!is.null(star_table)) {
      sel <- star_table[star_table$idx %in% selected_idx, , drop = FALSE]
      data.frame(
        segID = sel$segID,
        xcen  = sel$segstats_xcen,
        ycen  = sel$segstats_ycen,
        mag   = sel$segstats_mag,
        R50   = sel$segstats_R50,
        axrat = sel$fit_axrat,
        stringsAsFactors = FALSE)
    } else data.frame(),
    sky_pedestal_subtracted = as.numeric(sky_pedestal),
    pool_tier               = pool_tier,
    pool_tier_label          = pool_tier_label,
    psf_montage_before_profuse = psf_montage_drawn_early
  )
}

cat("[profuse_psf_v13.R] loaded.\n")
