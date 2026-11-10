#!/usr/bin/env Rscript
# Focused synthetic geometry proof, not an E1-E7/G5 gate or observed PSF fit.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('PSF rendering requires SLURM')
out <- Sys.getenv('BD_PSF_GEOMETRY_REPORT_DIR')
if (!nzchar(out) || !dir.exists(out)) stop('Missing fresh geometry report directory')
source('R/00_load_all.R')
root <- bd_find_root()
lib <- normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace('yaml',quietly=TRUE),requireNamespace('digest',quietly=TRUE),
          requireNamespace('jsonlite',quietly=TRUE),requireNamespace('ProFit',quietly=TRUE))
pin <- yaml::read_yaml(file.path(root,'config','package_pins.yaml'))$packages$ProFit
desc <- utils::packageDescription('ProFit',lib.loc=lib)
stopifnot(identical(desc$RemoteSha,pin$sha),
          package_version(desc$Version)==package_version(as.character(pin$version)),
          identical(normalizePath(getNamespaceInfo(asNamespace('ProFit'),'path')),
                    normalizePath(file.path(lib,'ProFit'))))
source(file.path(root,'R','psf_v14.R'))
hash <- function(path) digest::digest(file=path,algo='sha256')
dll <- getLoadedDLLs()[['ProFit']][['path']]
runtime <- list(R=R.version.string,ProFit_version=desc$Version,ProFit_sha=desc$RemoteSha,
  ProFit_path=find.package('ProFit'),ProFit_DLL=dll,ProFit_DLL_sha256=hash(dll),
  psf_v14_source_sha256=hash(file.path(root,'R','psf_v14.R')),
  diagnostic_source_sha256=hash(file.path(root,'tests','test_psf_centering.R')),
  git_HEAD=system2('git',c('rev-parse','HEAD'),stdout=TRUE),job=Sys.getenv('SLURM_JOB_ID'),
  scope='SYNTHETIC_GEOMETRY_ONLY_NO_E1_E7_OR_G5_CERTIFICATION')
print(runtime)
direct <- function(d,fwhm,con,q,ang) {
  ans <- ProFit::profitMakeModel(list(moffat=list(xcen=d/2,ycen=d/2,
    mag=0,fwhm=fwhm,con=con,axrat=q,ang=ang,box=0)),
    dim=as.integer(c(d,d)),psf=NULL,magzero=0)
  z <- if(is.list(ans)) ans$z else ans
  z/sum(z)
}
measure <- function(z,label,requested_dim) {
  stopifnot(is.matrix(z),all(is.finite(z)),all(z>=0),sum(z)>0)
  nr <- nrow(z);nc <- ncol(z);x <- seq_len(nr)-.5;y <- seq_len(nc)-.5
  cx <- sum(z*outer(x,rep(1,nc)))/sum(z)
  cy <- sum(z*outer(rep(1,nr),y))/sum(z)
  rotated <- z[nr:1,nc:1,drop=FALSE]
  row <- data.frame(case=label,requested_dim=requested_dim,nrow=nr,ncol=nc,
    sum=sum(z),physical_centroid_x=cx,physical_centroid_y=cy,
    physical_grid_centre_x=nr/2,physical_grid_centre_y=nc/2,
    delta_x=cx-nr/2,delta_y=cy-nc/2,
    normalized_180deg_asymmetry=max(abs(z-rotated))/max(z),
    object_sha256=digest::digest(z,algo='sha256'),stringsAsFactors=FALSE)
  row$centred <- max(abs(c(row$delta_x,row$delta_y)))<1e-6
  row$rotation_symmetric <- row$normalized_180deg_asymmetry<1e-6
  saveRDS(z,file.path(out,paste0(label,'_',requested_dim,'.rds')),version=2)
  print(row)
  row
}
rows <- list()
add <- function(z,label,d) rows[[length(rows)+1L]] <<- measure(z,label,d)
for(d in c(21L,31L,61L)) {
  for(shape in c('round','elliptical')) {
    q <- if(shape=='round') 1 else .67
    ang <- if(shape=='round') 0 else 37
    add(direct(d,3.2,2.7,q,ang),paste0('direct_physical_centre_',shape),d)
    add(.render_moffat(3.2,2.7,axrat=q,ang=ang,dim=c(d,d)),
        paste0('v14_render_',shape),d)
    fake <- rep(list(list(fit_par=list(fwhm=3.2,con=2.7,axrat=q,ang=ang))),3L)
    add(.build_median_param_psf(fake,max_dim=d,min_dim=d,
      median_segment_size_px=d,pixscale_arcsec=.3),paste0('v14_median_',shape),d)
  }
  # Broad enough to reach each requested cap, so all21/31/61 default sizes are tested.
  add(.build_default_psf(default_fwhm_px=6,default_con=2.1,
      psf_dim=c(d,d),pixscale_arcsec=.3),'v14_default_round',d)
}
metrics <- do.call(rbind,rows)
write.table(metrics,file.path(out,'centering_metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
controls <- grepl('^direct_',metrics$case)
bad <- !metrics$centred | !metrics$rotation_symmetric | abs(metrics$sum-1)>1e-12
report <- list(runtime=runtime,physical_coordinate_convention='Pixel i has centre i-0.5; grid centre is dim/2',
  controls_ok=all(!bad[controls]),current_v14_centred=all(!bad[!controls]),metrics=metrics,
  original_psf_code_modified=FALSE,
  scope='Geometry assertions only; original E1-E7/G5 criteria remain unspecified and unrun')
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,digits=NA),
           file.path(out,'centering_report.json'))
stopifnot(all(!bad[controls]))
if(any(bad[!controls])) stop('PSF_CENTERING_ASSERTIONS_FAILED; see recorded centroid and symmetry evidence')
cat('PSF_CENTERING_ASSERTIONS_SATISFIED_NO_SCIENTIFIC_GATE_PROMOTION\n')
