# Stage resume and cache provenance (Claude review of 148309f, 4 Nov 2026). Synthetic files only.
# 1. Resuming a saved stage must work (148309f: crashes, calls undefined .bd_stage_input_hash).
# 2. Observed mode: changing a verified index row (labels/coordinates) must invalidate the cache.
# 3. Observed mode: malformed source_files must fail closed, never fall back to cached hashes.
# 4. Observed mode: changed file bytes must invalidate the cache.
source(file.path('R','00_load_all.R')); bd_load_all()
config <- bd_config('mac')
out <- tempfile('resume_'); dir.create(out)
ctx <- bd_context('W999000001', config, run_id='resume_test', mode='synthetic', output_override=out)
invisible(bd_stage(ctx,'inputs'))
ctx2 <- bd_context('W999000001', config, run_id='resume_test', mode='synthetic', output_override=out)
r <- try(bd_stage(ctx2,'inputs'), silent=TRUE)
stopifnot('1. synthetic stage resume works' = !inherits(r,'try-error'))

td <- tempfile('idx_'); dir.create(td)
f <- function(name, bytes) { p <- file.path(td,name); writeBin(as.raw(bytes), p); p }
img <- f('img.fits',1:64); seg <- f('seg.fits',2:65); psf <- f('psf.rds',3:66)
row <- data.frame(uberID='W999000002', image_path=img, segmap_path=seg, image_hdu=1L, segmap_hdu=1L, psf_path=psf,
                  ra_deg=10, dec_deg=-30, segID=7L, psf_fwhm_px=2.2, detector_mask_policy='NONE_VERIFIED_FROM_SOURCE',
                  stringsAsFactors=FALSE)
index <- file.path(td,'verified_path_index.tsv'); utils::write.table(row, index, sep='\t', quote=FALSE, row.names=FALSE)
cfg <- config; cfg$paths$input_index <- index
octx <- bd_context('W999000002', cfg, run_id='prov_test', mode='observed', output_override=tempfile('o_'))
cached <- bd_index_row(cfg,'W999000002')
res <- list(result=list(observed=list(source_files=cached,
              source_hashes=list(image=cached$image_path_sha, segmap=cached$segmap_path_sha, psf=cached$psf_path_sha))))
h0 <- bd_stage_input_hash(res, octx)
stopifnot('baseline reproducible' = identical(h0, bd_stage_input_hash(res, octx)))
row2 <- row; row2$segID <- 99L; utils::write.table(row2, index, sep='\t', quote=FALSE, row.names=FALSE)
stopifnot('2. label change invalidates' = inherits(try(bd_stage_input_hash(res, octx), silent=TRUE),'try-error'))
utils::write.table(row, index, sep='\t', quote=FALSE, row.names=FALSE)
bad <- res; bad$result$observed$source_files <- NULL
stopifnot('3. malformed provenance fails closed' = inherits(try(bd_stage_input_hash(bad, octx), silent=TRUE),'try-error'))
bad_hash <- res
bad_hash$result$observed$source_hashes$psf <- NULL
stopifnot('missing PSF hash rejected' =
  inherits(try(bd_stage_input_hash(bad_hash,octx),silent=TRUE),'try-error'))
writeBin(as.raw(64:1), img)
stopifnot('4. changed bytes invalidate' = inherits(try(bd_stage_input_hash(res, octx), silent=TRUE),'try-error'))
cat('STAGE_RESUME_PROVENANCE_PASS\n')
