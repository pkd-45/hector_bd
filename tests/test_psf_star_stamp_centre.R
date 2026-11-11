#!/usr/bin/env Rscript
# Synthetic star-stamp coordinate convention for psf_v14 per-star Moffat fits.
# Supports the review patch hunks at .fit_one_star_gama_v11/v12 (stamp centre n/2).
# Not an E1-E7/G5 gate, no observed data and no change to suitability criteria.
#
# Facts checked on synthetic stars at known physical positions:
#  A1 ProFound segstats xcen/ycen use the physical convention (pixel i centre = i-0.5).
#  A2 With profuse_psf_v9.R sourced first (as scripts/build_pilot_psf_candidate.R does),
#     .refine_centroid returns 1-based index coordinates (pixel i centre = i).
#  A3 Without v9 the ProFound centre is used unchanged.
#  A4 Both star fitters recover the true centre whatever the initial guess.
#  A5 Production (v9) path: the star lies within 0.5 px of the stamp's physical centre n/2,
#     i.e. the patched initial guess and +/-3 px bounds are centred on the expected star.
# Reported, not asserted: .gama_suitability(_dr4) still tests |xcen-(n+1)/2| <= 1 px.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('PSF rendering requires SLURM')
out <- Sys.getenv('BD_PSF_GEOMETRY_REPORT_DIR')
if (!nzchar(out) || !dir.exists(out)) stop('Missing fresh geometry report directory')
source('R/00_load_all.R')
root <- bd_find_root()
lib <- normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace('yaml',quietly=TRUE),requireNamespace('digest',quietly=TRUE),
          requireNamespace('jsonlite',quietly=TRUE),requireNamespace('ProFit',quietly=TRUE),
          requireNamespace('ProFound',quietly=TRUE))
pins <- yaml::read_yaml(file.path(root,'config','package_pins.yaml'))$packages
for (package in c('ProFit','ProFound')) {
  desc <- utils::packageDescription(package,lib.loc=lib)
  stopifnot(identical(desc$RemoteSha,pins[[package]]$sha),
            identical(normalizePath(getNamespaceInfo(asNamespace(package),'path')),
                      normalizePath(file.path(lib,package))))
}
hash <- function(path) digest::digest(file=path,algo='sha256')
v9 <- Sys.getenv('BD_PSF_V9','/home/uqpdas/WAVES_BD_pipeline/profuse_psf_v9.R')
v14 <- file.path(root,'R','psf_v14.R')
stopifnot(file.exists(v9))
runtime <- list(R=R.version.string,job=Sys.getenv('SLURM_JOB_ID'),
  git_HEAD=system2('git',c('rev-parse','HEAD'),stdout=TRUE),
  ProFit_sha=utils::packageDescription('ProFit',lib.loc=lib)$RemoteSha,
  ProFound_sha=utils::packageDescription('ProFound',lib.loc=lib)$RemoteSha,
  psf_v14_source_sha256=hash(v14),psf_v9_path=normalizePath(v9),psf_v9_sha256=hash(v9),
  diagnostic_source_sha256=hash(file.path(root,'tests','test_psf_star_stamp_centre.R')),
  scope='SYNTHETIC_STAR_STAMP_CONVENTION_ONLY_NO_E1_E7_OR_G5_CERTIFICATION')
print(runtime)

environments <- list()
for (path_name in c('production_v9','without_v9')) {
  env <- new.env(parent=.GlobalEnv)
  if (path_name=='production_v9') sys.source(v9,envir=env)
  sys.source(v14,envir=env)
  stopifnot(identical(exists('.refine_centroid',envir=env,inherits=FALSE),path_name=='production_v9'))
  environments[[path_name]] <- env
}

dim_img <- c(81L,81L); label <- 7L; base <- 40.5   # pixel index 41 spans [40,41]
offsets <- data.frame(dx=c(-0.45,-0.25,0,0.25,0.45,-0.45,0.45),
                      dy=c(-0.45,-0.25,0,0.25,0.45,0.45,-0.45))
rows <- list()
set.seed(20261111L)
for (k in seq_len(nrow(offsets))) {
  x0 <- base+offsets$dx[k]; y0 <- base+offsets$dy[k]
  model <- ProFit::profitMakeModel(list(moffat=list(xcen=x0,ycen=y0,mag=-10,fwhm=3.2,
    con=2.7,axrat=1,ang=0,box=0)),dim=dim_img,magzero=0,psf=NULL)
  z <- if (is.list(model)) model$z else model
  image <- z+matrix(stats::rnorm(prod(dim_img)),dim_img[1],dim_img[2])
  px <- seq_len(dim_img[1])-0.5; py <- seq_len(dim_img[2])-0.5
  segim <- matrix(0L,dim_img[1],dim_img[2])
  segim[sqrt(outer((px-x0)^2,(py-y0)^2,'+'))<=6] <- label
  stats_table <- ProFound::profoundSegimStats(image=image,segim=segim,sky=0,skyRMS=1,
    magzero=0,pixscale=1)
  seg_row <- stats_table[stats_table$segID==label,,drop=FALSE]
  stopifnot(nrow(seg_row)==1L)
  for (path_name in names(environments)) {
    env <- environments[[path_name]]
    for (fitter in c('v11_BFGS','v12_LBFGSB')) {
      fn <- get(if (fitter=='v11_BFGS') '.fit_one_star_gama_v11' else '.fit_one_star_gama_v12',
                envir=env,inherits=FALSE)
      f <- fn(image,segim,seg_row=seg_row,half_size=12L,magzero=0)
      if (!isTRUE(f$ok)) stop(sprintf('Star fit failed: %s %s case %d: %s',
                                      path_name,fitter,k,paste(f$reason,collapse=' ')))
      n <- f$stamp_size
      tx <- x0-(f$cutout_anchor[['x_lo']]-1); ty <- y0-(f$cutout_anchor[['y_lo']]-1)
      suit <- get('.gama_suitability_dr4',envir=env,inherits=FALSE)(list(f),
        do_2sig_clip=FALSE,verbose=FALSE)
      rows[[length(rows)+1L]] <- data.frame(path=path_name,fitter=fitter,case=k,
        true_dx_from_pixel_centre=offsets$dx[k],true_dy_from_pixel_centre=offsets$dy[k],
        profound_xcen_minus_true=seg_row$xcen-x0,profound_ycen_minus_true=seg_row$ycen-y0,
        stamp_xc_minus_true=f$segstats_xcen-x0,stamp_yc_minus_true=f$segstats_ycen-y0,
        stamp_size=n,true_stamp_x=tx,true_stamp_y=ty,
        fit_x_error=f$fit_par$xcen-tx,fit_y_error=f$fit_par$ycen-ty,
        true_minus_patched_centre_x=tx-n/2,true_minus_patched_centre_y=ty-n/2,
        true_minus_unpatched_centre_x=tx-(n+1)/2,true_minus_unpatched_centre_y=ty-(n+1)/2,
        suitability_offset_x=f$fit_par$xcen-(n+1)/2,suitability_offset_y=f$fit_par$ycen-(n+1)/2,
        suitability_centre_ok=abs(f$fit_par$xcen-(n+1)/2)<=1&&abs(f$fit_par$ycen-(n+1)/2)<=1,
        suitability_dr4_all_criteria=isTRUE(suit[1]),reduced_chisq=f$reduced_chisq,
        fit_mag=f$fit_par$mag,stringsAsFactors=FALSE)
    }
  }
}
metrics <- do.call(rbind,rows)
write.table(metrics,file.path(out,'star_stamp_centre_metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(metrics[,c('path','fitter','case','stamp_xc_minus_true','fit_x_error','fit_y_error',
                 'true_minus_patched_centre_x','true_minus_unpatched_centre_x',
                 'suitability_offset_x','suitability_offset_y','suitability_centre_ok')])
prod <- metrics$path=='production_v9'
checks <- list(
  A1_profound_physical_convention=all(abs(c(metrics$profound_xcen_minus_true,
    metrics$profound_ycen_minus_true))<0.05),
  A2_v9_refine_index_convention=all(abs(c(metrics$stamp_xc_minus_true[prod],
    metrics$stamp_yc_minus_true[prod])-0.5)<0.1),
  A3_without_v9_uses_profound_centre=all(abs(c(metrics$stamp_xc_minus_true[!prod],
    metrics$stamp_yc_minus_true[!prod]))<0.05),
  A4_fit_recovers_true_centre=all(abs(c(metrics$fit_x_error,metrics$fit_y_error))<0.05),
  A5_production_star_within_half_pixel_of_patched_centre=all(abs(c(
    metrics$true_minus_patched_centre_x[prod],metrics$true_minus_patched_centre_y[prod]))<=0.6))
offsets_prod <- c(metrics$suitability_offset_x[prod],metrics$suitability_offset_y[prod])
report <- list(runtime=runtime,checks=checks,
  production_expected_star_centre_in_stamp='n/2 (+/-0.5): v9 .refine_centroid returns index coordinates',
  without_v9_expected_star_centre_in_stamp='(n+1)/2 (+/-0.5): ProFound physical centre is used directly',
  suitability_window_centre_in_code='(stamp_size+1)/2 (R/psf_v14.R .gama_suitability and .gama_suitability_dr4)',
  production_suitability_offset_range=range(offsets_prod),
  production_suitability_centre_failures=sum(!metrics$suitability_centre_ok[prod]),
  suitability_window_note=paste('In the production path the 1 px window spans -0.5..+1.5 px about the',
    'expected star centre. Reported for scientific review; the criterion is not changed here.'),
  metrics=metrics)
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,digits=NA),
           file.path(out,'star_stamp_centre_report.json'))
print(checks)
if (!all(unlist(checks))) stop('PSF_STAR_STAMP_CONVENTION_ASSERTIONS_FAILED')
cat('PSF_STAR_STAMP_CONVENTION_ASSERTIONS_SATISFIED_NO_SCIENTIFIC_GATE_PROMOTION\n')
