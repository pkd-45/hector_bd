# D7: preflight DESIGN ONLY, no optimiser, Highlander, MCMC or real fit call.
# User must approve after Step0 evidence and the 20-galaxy synthetic/pilot gates.
# D7 real-data execution is deliberately absent from this review candidate.
bd_d7_preflight <- function(built,originals,models,tofit,tolog,intervals,
                            constraints,magzero,finesample,gate_evidence) {
  required <- c('G0_STEP0','G2_SUPPORT','G3_HASHES','G4_F1_F9','G5_PSF','G6_MODELS')
  if (!is.list(gate_evidence) || !setequal(names(gate_evidence),required) ||
      !all(vapply(gate_evidence,isTRUE,logical(1))))
    stop('G0-G6 not independently approved; D7 cannot proceed')
  if (!is.list(originals) || !all(c('L','image','sigma','detector_mask') %in% names(originals)))
    stop('Independently reloaded original arrays required for G2/G3')
  bd_assert_inputs_unchanged(built,originals$L,originals$image,
                             originals$sigma,originals$detector_mask)
  if (!is.numeric(magzero) || length(magzero)!=1L || !is.finite(magzero))
    stop('Explicit verified magzero required')
  if (!(finesample %in% c(1L,3L,5L)) || length(finesample)!=1L)
    stop('Explicit finesample=1,3,5 required')
  families <- c('S','P','D_FIXED4','BD_FREE_N')
  for (n in families)
    if (is.null(models[[n]]) || is.null(tofit[[n]]) ||
        is.null(tolog[[n]]) || is.null(intervals[[n]]) ||
        is.null(constraints[[n]]))
      stop('Model family lacks frozen definition, bounds or constraints: ',n)
  list(status='D7_PREFLIGHT_ONLY_NO_FIT',families=families,
       original_support_hash=built$contract$hashes$support,
       original_segmap_hash=built$contract$hashes$L,
       input_hashes=built$hashes,magzero=magzero,finesample=finesample)
}
