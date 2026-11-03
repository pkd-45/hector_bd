# D2: original WAVES target-pixel contract. Pure functions, no FITS IO or fitting.
# All R positions are (row, column), one-based, in the ORIGINAL science-image grid.
.bd_assert <- function(ok, why) if (!isTRUE(ok)) stop(why, call. = FALSE)
.bd_hash <- function(x) {
  if (!requireNamespace('digest', quietly = TRUE)) stop('Requires R package digest')
  digest::digest(x, algo = 'sha256', serialize = TRUE)
}
.bd_required_step0 <- c('segmentation_product', 'mapping_evidence', 'science_band',
  'image_units', 'sky_state', 'science_image_path', 'segmap_path',
  'science_hdu', 'segmap_hdu', 'sigma_source', 'gain_source', 'wcs_verified')
validate_step0 <- function(step0) {
  .bd_assert(is.list(step0), 'Step0 must be a named list, not a default')
  missing <- setdiff(.bd_required_step0, names(step0))
  .bd_assert(!length(missing), paste('Missing Step0 fields:', paste(missing, collapse = ',')))
  for (n in .bd_required_step0) {
    .bd_assert(length(step0[[n]]) == 1L && !is.na(step0[[n]]) &&
      nzchar(trimws(as.character(step0[[n]]))), paste('Unresolved Step0:', n))
  }
  .bd_assert(step0$segmentation_product %in% c('detection','dilated_total','regrouped','as_delivered_unclassified'),
    'Declare delivered map product or as_delivered_unclassified; original segmap path and header must be verified')
  .bd_assert(identical(step0$wcs_verified, TRUE), 'Original image/segmap grid and WCS need independent verification')
  invisible(TRUE)
}
segment_contract <- function(L, A_g, image, sigma, detector_mask, step0) {
  validate_step0(step0)
  .bd_assert(is.matrix(L) && (is.integer(L) || is.numeric(L)), 'L must be original label matrix')
  .bd_assert(all(is.finite(L)) && all(L == floor(L)) && all(L >= 0), 'L must have finite nonnegative integer labels')
  .bd_assert(is.numeric(A_g) && length(A_g) > 0L && all(is.finite(A_g)) &&
    all(A_g > 0) && all(A_g == floor(A_g)) && !anyDuplicated(A_g), 'A_g must be documented distinct positive integer labels')
  .bd_assert(all(A_g %in% L), 'Mapped label missing in original segmentation')
  .bd_assert(is.matrix(image) && is.matrix(sigma) && is.matrix(detector_mask), 'Image, sigma and detector mask must be matrices')
  .bd_assert(identical(dim(image),dim(L)) && identical(dim(sigma),dim(L)) &&
    identical(dim(detector_mask),dim(L)), 'Original image/grid mismatch; no resampling permitted')
  .bd_assert(is.numeric(image) && is.numeric(sigma), 'Image and sigma must be numeric')
  .bd_assert(is.logical(detector_mask) && !anyNA(detector_mask),
    'Detector mask must be logical: TRUE=bad, with independently documented provenance')
  T_g <- L %in% A_g
  dim(T_g) <- dim(L)
  bad_image <- !is.finite(image)
  bad_sigma <- !is.finite(sigma) | sigma <= 0
  bad_detector <- detector_mask
  Q <- !(bad_image | bad_sigma | bad_detector)
  F_g <- T_g & Q
  .bd_assert(any(T_g), 'No target pixels; catalogue-to-label mapping unresolved')
  .bd_assert(any(F_g), 'All target pixels unusable; do not fit')
  # No connected-component filtering, clipping, morphology, or label edits.
  .bd_assert(identical(F_g,T_g & Q), 'G2: original-grid support identity failed')
  dropped <- which(T_g & !Q, arr.ind = TRUE)
  ledger <- data.frame(row_1based = dropped[,1], col_1based = dropped[,2],
    row_y0 = dropped[,1] - 1L, col_x0 = dropped[,2] - 1L,
    cause = ifelse(bad_image[T_g & !Q], 'NONFINITE_IMAGE',
      ifelse(bad_sigma[T_g & !Q], 'INVALID_SIGMA', 'DETECTOR_MASK')),
    bad_image = bad_image[T_g & !Q], bad_sigma = bad_sigma[T_g & !Q],
    detector_mask = bad_detector[T_g & !Q], stringsAsFactors = FALSE)
  .bd_assert(nrow(ledger) == sum(T_g) - sum(F_g), 'Ledger incomplete')
  coords <- which(F_g, arr.ind = TRUE)
  list(T_g = T_g, Q = Q, F_g = F_g, ledger = ledger,
    n_original_target = sum(T_g), n_usable = sum(F_g),
    hashes = list(L = .bd_hash(L), target = .bd_hash(T_g), Q = .bd_hash(Q),
      support = .bd_hash(F_g), support_coords = .bd_hash(coords),
      image = .bd_hash(image), sigma = .bd_hash(sigma), detector_mask = .bd_hash(detector_mask)),
    labels = as.integer(A_g), step0 = step0, grid = dim(L))
}
assert_segment_g2 <- function(contract, L, image, sigma, detector_mask, original_region) {
  .bd_assert(identical(dim(L), contract$grid), 'G2 original dimension changed')
  .bd_assert(identical(original_region, contract$F_g), 'G2 support differs from original WAVES label and usability')
  .bd_assert(identical(.bd_hash(L), contract$hashes$L), 'Original segmentation mutated')
  .bd_assert(identical(.bd_hash(image), contract$hashes$image), 'Science image mutated')
  .bd_assert(identical(.bd_hash(sigma), contract$hashes$sigma), 'Sigma mutated')
  .bd_assert(identical(.bd_hash(detector_mask), contract$hashes$detector_mask), 'Detector mask mutated')
  .bd_assert(identical(.bd_hash(original_region), contract$hashes$support), 'Support mutated')
  invisible(TRUE)
}

# A1: resolve label(s) by source catalogue RA/Dec using a VERIFIED converter
# for the ORIGINAL segmap WCS. The converter must return one-based R row/col
# integer indices; this function never uses nearest-labelled-pixel fallback.
# Pixel centres are checked in a 3x3 label neighbourhood (including sky).
bd_resolve_original_labels <- function(L, ra_deg, dec_deg, segmap_wcs_to_index,
   catalogue_seg_id=NULL, catalogue_group_id=NULL,
   group_members=NULL, verified_group_id=NULL,
   catalogue_mag=NULL, segment_mag=NULL) {
  .bd_assert(is.matrix(L) && all(is.finite(L)) && all(L==floor(L)) && all(L>=0),
    'A1 requires original nonnegative integer segmentation')
  .bd_assert(is.finite(ra_deg) && length(ra_deg)==1L &&
    is.finite(dec_deg) && length(dec_deg)==1L, 'A1 catalogue coordinates missing')
  .bd_assert(is.function(segmap_wcs_to_index) &&
    identical(attr(segmap_wcs_to_index,'verified_original_segmap_wcs'),TRUE),
    'A1 requires a verified original segmap-WCS converter')
  pos <- segmap_wcs_to_index(ra_deg,dec_deg)
  .bd_assert(is.list(pos) && identical(pos$frame,'ORIGINAL_SEGMAP') &&
    identical(pos$convention,'R_1BASED_ROW_COL'),
    'A1 converter must explicitly return original segmap R indices')
  rr <- pos$row; cc <- pos$col
  .bd_assert(length(rr)==1L && length(cc)==1L && is.finite(rr) && is.finite(cc) &&
    rr==floor(rr) && cc==floor(cc), 'A1 WCS pixel indices invalid')
  rr <- as.integer(rr); cc <- as.integer(cc)
  .bd_assert(rr>=1L && rr<=nrow(L) && cc>=1L && cc<=ncol(L),
    'A1 WCS location outside original map')
  label <- as.integer(L[rr,cc])
  .bd_assert(label>0L, 'A1 BLOCKED: WCS catalogue position is sky (L=0); no nearest-label fallback')
  .bd_assert(rr>1L && rr<nrow(L) && cc>1L && cc<ncol(L),
    'A1 BLOCKED: 3x3 boundary check unavailable at image edge')
  boundary_flag <- !all(L[(rr-1L):(rr+1L),(cc-1L):(cc+1L)]==label)
  # Amendment A1: boundary is FLAGGED, but the exact centre pixel label is used.
  # Never replace it with a neighbour or enlarge the membership.
  if (!is.null(catalogue_seg_id)) {
    .bd_assert(length(catalogue_seg_id)==1L && is.finite(catalogue_seg_id) &&
      catalogue_seg_id==label,'A1 BLOCKED: catalogue segID contradicts original label')
  }
  if (!is.null(catalogue_group_id)) {
    .bd_assert(!is.null(verified_group_id) &&
      identical(as.character(verified_group_id),as.character(catalogue_group_id)),
      'A1 BLOCKED: catalogue groupID not verified against group-membership table')
    .bd_assert(!is.null(group_members) && length(group_members)>=1L &&
      is.numeric(group_members) && all(is.finite(group_members)) &&
      all(group_members>0) && all(group_members==floor(group_members)) &&
      !anyDuplicated(group_members) && label %in% group_members &&
      all(group_members%in%L),
      'A1 BLOCKED: groupID members do not contain valid original labels')
    A_g <- sort(as.integer(group_members))
  } else {
    .bd_assert(is.null(group_members),
      'A1 BLOCKED: group members provided without catalogue groupID')
    A_g <- label
  }
  px <- which(matrix(L %in% A_g, nrow(L), ncol(L)), arr.ind=TRUE)
  centroid_offset <- sqrt((mean(px[,1])-rr)^2+(mean(px[,2])-cc)^2)
  mag_diff <- if(!is.null(catalogue_mag) && !is.null(segment_mag) &&
     length(catalogue_mag)==1L && length(segment_mag)==1L &&
     is.finite(catalogue_mag) && is.finite(segment_mag))
     as.numeric(segment_mag-catalogue_mag) else NA_real_
  list(A_g=A_g, row=rr,col=cc,central_label=label,
    boundary_check=if (boundary_flag) 'FLAGGED_3x3_LABEL_BOUNDARY' else '3x3_UNIFORM_LABEL',
    boundary_flag=boundary_flag,
    centroid_offset_px=centroid_offset,mag_delta=mag_diff,
    provenance='VERIFIED_SEGMAP_WCS_POSITION_AND_AVAILABLE_CATALOGUE_CHECKS')
}
