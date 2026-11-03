#!/usr/bin/env python3
"""Deterministic patch of the byte-verified historical v13 source; no fitting."""
from pathlib import Path
import hashlib,sys
old=Path(sys.argv[1]); new=Path(sys.argv[2]); s=old.read_text()
EXPECTED='1de310c015d027c25b90f0898cfc7011e0ab8ef3e30d889835db2d3dd569e9f0'
assert hashlib.sha256(old.read_bytes()).hexdigest()==EXPECTED, 'V13 SOURCE HASH MISMATCH - refuse to patch'
assert s.count('stats::mad(v[keep]) * 1.4826')==1
assert s.count('stats::mad(sky_pix2, na.rm = TRUE) * 1.4826')==2
assert s.count('stats::mad(as.numeric(data_stamp), na.rm = TRUE) * 1.4826')==2
s=s.replace('stats::mad(v[keep]) * 1.4826','stats::mad(v[keep])')
s=s.replace('stats::mad(sky_pix2, na.rm = TRUE) * 1.4826','stats::mad(sky_pix2, na.rm = TRUE)')
s=s.replace('stats::mad(as.numeric(data_stamp), na.rm = TRUE) * 1.4826','stats::mad(as.numeric(data_stamp), na.rm = TRUE)')
# Two DR4 and legacy PA gates; axial angles modulo 180 degrees, not linear on [-90,90).
assert s.count('.iter_2sig_clip(angs,  in_log = FALSE)')==2
s=s.replace('.iter_2sig_clip(angs,  in_log = FALSE)','.v14_axial_pa_clip(angs)')
anchor='.iter_2sig_clip <- function(x, in_log = FALSE, sigma = 2, max_iter = 5) {'
helper='''# v14: PA is axial (period 180 degrees), so +/-89 degrees are adjacent.
.v14_axial_pa_clip <- function(x, sigma=2, max_iter=5) {
  ok <- is.finite(x)
  if (sum(ok)<3L) return(ok)
  ang <- x[ok]*pi/90  # doubled PA, radians
  keep <- rep(TRUE,length(ang))
  for (it in seq_len(max_iter)) {
    center <- atan2(stats::median(sin(ang[keep])),
                    stats::median(cos(ang[keep])))
    residual <- atan2(sin(ang-center), cos(ang-center))/2
    spread <- stats::mad(residual[keep])
    if (!is.finite(spread) || spread<=0) break
    nxt <- abs(residual) < sigma*spread
    if (identical(nxt,keep)) break
    keep <- nxt
    if (sum(keep)<3L) break
  }
  out <- ok;out[which(ok)] <- keep;out
}

'''
assert s.count(anchor)==1
s=s.replace(anchor,helper+anchor)
# Track flux loss against analytic infinite Moffat model and square/corner crop.
# 720-direction quadrature computes actual elliptical Moffat mass INSIDE an
# inscribed circle; since zero_corners removes the exterior, this is a
# physically meaningful expectation, not an exact pixel-sampled flux ratio.
anchor='.build_median_param_psf <- function(selected_fits,'
ee_helper='''# v14: approximate encircled fraction for an elliptical Moffat in a circular
# aperture; midpoint angular quadrature of the analytic radial integral.
# This must be labelled ANALYTIC_APPROX, never mistaken for measured star flux.
.v14_moffat_circle_ee <- function(fwhm,beta,q,diameter,nang=720L) {
  if (!all(is.finite(c(fwhm,beta,q,diameter))) || fwhm<=0 || beta<=1 ||
      q<=0 || q>1 || diameter<1) return(NA_real_)
  alpha <- (fwhm/2)/sqrt(2^(1/beta)-1)
  theta <- 2*pi*((seq_len(nang)-0.5)/nang)
  v <- cos(theta)^2 + (sin(theta)/q)^2
  tail <- mean((1/v)*(1 + v*(diameter/2/alpha)^2)^(1-beta))/q
  pmax(0,pmin(1,1-tail))
}
.v14_annotate_ee <- function(psf,fwhm,beta,q,min_fraction=0.99) {
  if (is.null(psf)) return(NULL)
  fraction <- .v14_moffat_circle_ee(fwhm,beta,q,min(dim(psf)))
  attr(psf,'ee_infinite_moffat_approx') <- fraction
  attr(psf,'ee_method') <- 'ELLIPTICAL_MOFFAT_ANALYTIC_CIRCULAR_APERTURE_720_ANGLES_APPROX'
  attr(psf,'ee_lost_fraction_approx') <- if (is.finite(fraction)) 1-fraction else NA_real_
  attr(psf,'ee_below_99pct') <- !is.finite(fraction) || fraction<min_fraction
  attr(psf,'ee_decision') <- if (!is.finite(fraction)) 'UNRESOLVED'
      else if (fraction<min_fraction) 'LOSS_RECORDED_NOT_SILENT'
      else '>=99_PERCENT_APPROX'
  psf
}

'''
assert s.count(anchor)==1;s=s.replace(anchor,ee_helper+anchor)
# Both builders use corner mask; record the analytic fraction explicitly.
old='''  psf <- .zero_corners(psf)

  attr(psf, "fwhm")'''
assert s.count(old)==1
s=s.replace(old,'''  psf <- .zero_corners(psf)
  psf <- .v14_annotate_ee(psf, fwhm_med, con_med, axr_med)

  attr(psf, "fwhm")''',1)
s=s.replace('  psf <- .zero_corners(psf)\n  attr(psf, "fwhm")',
'''  psf <- .zero_corners(psf)
  psf <- .v14_annotate_ee(psf, default_fwhm_px, default_con, 1)
  attr(psf, "fwhm")''',1)
# Mark non-analytic selected kernels as unresolved rather than silently
# recording a null EE number. Insert AFTER kind inference: changing attrs
# earlier would break v13's identical() dispatch.
kind_end='''                      else                                                              "default"
  if (verbose)'''
assert s.count(kind_end)==2
s=s.replace(kind_end,'''                      else                                                              "default"
  if (is.null(attr(recommended, 'ee_decision'))) {
    attr(recommended, 'ee_infinite_moffat_approx') <- NA_real_
    attr(recommended, 'ee_lost_fraction_approx') <- NA_real_
    attr(recommended, 'ee_method') <- 'NONANALYTIC_PSF_EE_UNMEASURED'
    attr(recommended, 'ee_decision') <- 'BLOCKED_MEASURE_PSF_EE_BEFORE_PILOT'
  }
  if (verbose)''')
# add explicit per-run metadata to top-level result; metadata is also retained
# on the PSF array if consumers save only the kernel.
needle='''    output_kind            = output_kind_used,'''
assert s.count(needle)==2
s=s.replace(needle,needle+'''
    ee_infinite_moffat_approx = attr(recommended, 'ee_infinite_moffat_approx'),
    ee_lost_fraction_approx = attr(recommended, 'ee_lost_fraction_approx'),
    ee_method = attr(recommended, 'ee_method'),
    ee_decision = attr(recommended, 'ee_decision'),''')
# Do not overwrite global .render_moffat normalisation behavior (historical
# ProFit code). The loss is no longer SILENT: all analytic variants get attrs.
assert 'stats::mad(v[keep]) * 1.4826' not in s
assert s.count('.v14_axial_pa_clip(angs)')==2
s = '# v14 generated from byte-verified historical v13. Synthetic-only review candidate.\n' + s
new.write_text(s)
print('v14 patched:',new,hashlib.sha256(new.read_bytes()).hexdigest())
