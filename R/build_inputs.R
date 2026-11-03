# D3: deterministic, zero-fit input builder. Requires segment_contract.R sourced first.
# Source data are already on one VERIFIED original pixel grid; no regridding occurs.
bd_build_inputs <- function(L, A_g, image, sigma, detector_mask, psf, fwhm_px, step0) {
  contract <- segment_contract(L, A_g, image, sigma, detector_mask, step0)
  .bd_assert(is.matrix(psf) && is.numeric(psf) && nrow(psf) == ncol(psf) &&
    nrow(psf) %% 2L == 1L && all(is.finite(psf)) && all(psf >= 0) &&
    sum(psf) > 0, 'PSF must be a finite nonnegative, square odd-sized matrix')
  .bd_assert(length(fwhm_px) == 1L && is.finite(fwhm_px) && fwhm_px > 0,
    'FWHM must come from declared PSF provenance in original pixels')
  # Pre-normalise ONCE before ProFit (older libprofit mutates caller PSF).
  psf_normalised <- psf/sum(psf)
  .bd_assert(abs(sum(psf_normalised)-1) < 1e-12, 'PSF normalisation failed')
  pos <- which(contract$T_g, arr.ind = TRUE)
  rr <- range(pos[,1]); cc <- range(pos[,2])
  pad <- as.integer(max(ceiling(nrow(psf_normalised)/2), ceiling(3*fwhm_px)))
  y <- seq.int(max(1L,rr[1]-pad),min(nrow(L),rr[2]+pad))
  x <- seq.int(max(1L,cc[1]-pad),min(ncol(L),cc[2]+pad))
  # The stamp padding is VISUAL transparency, not a convolution correctness fix.
  region <- contract$F_g[y,x,drop=FALSE]
  # Original-mask outside region is not modified: no autoclip, island removal,
  # cutbox, re-detection, relabelling, dilation or smoothing.
  input <- list(image=image[y,x,drop=FALSE], sigma=sigma[y,x,drop=FALSE],
    region=region, segim=L[y,x,drop=FALSE], mask=detector_mask[y,x,drop=FALSE],
    psf=psf_normalised)
  g2_full <- matrix(FALSE,nrow(L),ncol(L));g2_full[y,x] <- input$region
  .bd_assert(identical(g2_full,contract$F_g), 'G2 failed: cropped region did not reconstruct original support')
  assert_segment_g2(contract,L,image,sigma,detector_mask,g2_full)
  hashes <- lapply(input,.bd_hash)
  list(input=input,contract=contract,hashes=hashes,
    coords=list(row_1based=range(y),col_1based=range(x),
      row0=as.integer(min(y)-1L),col0=as.integer(min(x)-1L),
      pad_requested=pad,source_edge=any(c(rr[1]==1L,cc[1]==1L,
       rr[2]==nrow(L),cc[2]==ncol(L))),
      pad_truncated=any(c(min(y)==1L,max(y)==nrow(L),min(x)==1L,max(x)==ncol(L)))),
    sigma_mode='PROVIDED_VERIFIED_NO_POISSON_TERM_IN_INPUT_BUILDER',
    sky_mode='IMAGE_AS_VERIFIED_SKY_STATE_NO_AUTOMATIC_SUBTRACTION')
}
bd_assert_inputs_unchanged <- function(built, L, image, sigma, detector_mask) {
  assert_segment_g2(built$contract,L,image,sigma,detector_mask,
                    built$contract$F_g)
  for (n in names(built$input))
    .bd_assert(identical(.bd_hash(built$input[[n]]),built$hashes[[n]]),
      paste('G3 input mutated:',n))
  invisible(TRUE)
}
