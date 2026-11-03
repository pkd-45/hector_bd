# Selection/effective PSF comparison AFTER v13 and v14 were run on the same
# EXACT raw image, star catalogue, source segmap and explicit settings.
# This is deliberately a reporting helper, NOT a call to either builder.
.bd_psf_ids <- function(x) {
  stopifnot(is.list(x),!is.null(x$selected_idx),!is.null(x$fits))
  vapply(x$fits[x$selected_idx],function(z) as.character(z$segID),character(1))
}
compare_psf_selection <- function(v13_result,v14_result) {
  a <- .bd_psf_ids(v13_result); b <- .bd_psf_ids(v14_result)
  list(v13_star_ids=a,v14_star_ids=b,removed_from_v13=setdiff(a,b),
       added_by_v14=setdiff(b,a),identical_star_ids=identical(a,b),
       kind_v13=v13_result$output_kind,kind_v14=v14_result$output_kind,
       ee_v14=v14_result$ee_infinite_moffat_approx,
       ee_loss_v14=v14_result$ee_lost_fraction_approx,
       note='NOT A REAL OBSERVATION unless exact input hashes match')
}
