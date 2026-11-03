# Phase C diagnostic-only background mask, separate from the immutable L.
# A sky exclusion mask may be dilated; the galaxy segmentation/fit region never is.
bd_sky_exclusion <- function(L, radius_px) {
  stopifnot(is.matrix(L),length(radius_px)==1L,is.finite(radius_px),
            radius_px>=0, radius_px==as.integer(radius_px))
  n <- nrow(L);m <- ncol(L)
  excluded <- L!=0L
  original <- excluded
  radius <- as.integer(radius_px)
  if(radius>0L) {
    # O(radius * pixels), exact disc dilation, without an R morphology package.
    dist_horizontal <- matrix(Inf,n,m)
    for(j in seq_len(m)) {
      dist_horizontal[,j] <- ifelse(original[,j],0,
        if(j==1L) Inf else dist_horizontal[,j-1L]+1)
    }
    for(j in seq.int(m,1L)) {
      if(j<m) dist_horizontal[,j] <- pmin(dist_horizontal[,j],
                                           dist_horizontal[,j+1L]+1)
    }
    for(dx in seq.int(-radius,radius)) {
      lo <- max(1L,1L-dx); hi <- min(n,n-dx)
      if(lo>hi) next
      src <- seq.int(lo,hi)
      dst <- src+dx
      width <- floor(sqrt(radius*radius-dx*dx))
      excluded[dst,] <- excluded[dst,,drop=FALSE] |
                         (dist_horizontal[src,,drop=FALSE]<=width)
    }
  }
  excluded
}
bd_sky_noise <- function(image,L,effective_gain,sky_policy='fixed_measured',
                         sky_mask_radius_px=10L) {
  stopifnot(identical(dim(image),dim(L)),is.finite(effective_gain),effective_gain>0)
  if(!identical(sky_policy,'fixed_measured'))
    stop('Only measured fixed sky supported in primary pipeline')
  # Deliberately hold a separate mask; L is not modified.
  before <- .bd_hash(L)
  sky_exclusion <- bd_sky_exclusion(L,sky_mask_radius_px)
  if(!identical(before,.bd_hash(L))) stop('Sky computation mutated the original segmap')
  sky_pixels <- !sky_exclusion & L==0L & is.finite(image)
  b <- image[sky_pixels]
  if(length(b)<100L) stop('Too few clean sky pixels after mask dilation')
  med <- stats::median(b)
  rms <- 1.4826*stats::median(abs(b-med))
  if(!is.finite(rms) || rms<=0) stop('Sky RMS invalid')
  raw <- image[L==0L & is.finite(image)]
  rawmed <- stats::median(raw)
  sigma <- matrix(rms,nrow(image),ncol(image))
  list(residual_sky=med,skyRMS=rms,gain=effective_gain,sigma=sigma,
       noise='MODEL_DEPENDENT_SIGMA_IN_FIT; STATIC_FLOOR_FOR_G2',
       background_n=length(b),background_kind='OUTSIDE_SEPARATE_DILATED_SKY_EXCLUSION',
       sky_mask_radius_px=sky_mask_radius_px,
       sky_mask_sha256=.bd_hash(sky_exclusion),
       original_L_sha256=before,raw_L0_median=rawmed,
       median_shift_sigma=(rawmed-med)/rms,
       n_excluded=sum(sky_exclusion))
}
bd_model_sigma <- function(model,skyRMS,gain) {
  sqrt(skyRMS^2+pmax(model,0)/gain)
}
