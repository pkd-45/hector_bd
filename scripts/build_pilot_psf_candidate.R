#!/usr/bin/env Rscript
# One-target PSF-star measurements. NO observed galaxy fit or E1-E7/G5 decision.
if(!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('PSF-star fitting requires SLURM')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop('usage: build_pilot_psf_candidate.R validated_input_dir NEW_output_dir')
source('R/00_load_all.R')
root <- bd_find_root()
lib <- normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace('yaml',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE),
          requireNamespace('digest',quietly=TRUE))
source('R/config_io.R');source('R/io_fits.R')
input <- normalizePath(args[1],mustWork=TRUE)
out <- normalizePath(args[2],mustWork=FALSE)
if(!startsWith(input,paste0(root,'/outputs/pilot_inputs_')) ||
   !startsWith(out,paste0(root,'/outputs/')) || dir.exists(out) || file.exists(out))
  stop('Require validated repository pilot_inputs directory and NEW outputs directory')
manifest_path <- file.path(input,'copy_manifest.json')
verification_path <- file.path(input,'Rfits_verification.json')
after_path <- file.path(input,'originals_after.json')
stopifnot(file.exists(manifest_path),file.exists(verification_path),file.exists(after_path))
manifest <- jsonlite::fromJSON(manifest_path,simplifyVector=FALSE)
verification <- jsonlite::fromJSON(verification_path,simplifyVector=FALSE)
after <- jsonlite::fromJSON(after_path,simplifyVector=FALSE)
stopifnot(identical(manifest$uberID,'W11690092907692'),isTRUE(after$originals_unchanged),
          manifest$segID==4211L,manifest$ra_deg==11.091131080899515,
          manifest$dec_deg== -31.1410252394602)
lock <- yaml::read_yaml(file.path(root,'config','package_pins.yaml'))
package_manifest <- lapply(names(lock$packages),function(package) {
  pin <- lock$packages[[package]]
  desc <- utils::packageDescription(package,lib.loc=lib)
  stopifnot(identical(desc$RemoteSha,pin$sha))
  if(!is.null(pin$version)) stopifnot(package_version(desc$Version)==package_version(as.character(pin$version)))
  loadNamespace(package,lib.loc=lib)
  path <- normalizePath(getNamespaceInfo(asNamespace(package),'path'),mustWork=TRUE)
  stopifnot(identical(path,normalizePath(file.path(lib,package),mustWork=TRUE)))
  list(package=package,version=desc$Version,source_sha=desc$RemoteSha,path=path)
})
if(as.character(getRversion())!='4.3.3') stop('Expected verified R4.3.3 runtime')
hash <- function(path) digest::digest(file=path,algo='sha256')
object_hash <- function(x) digest::digest(x,algo='sha256')
numeric_hash <- function(x) digest::digest(writeBin(as.double(x),raw(),size=8L,
  endian='little'),algo='sha256',serialize=FALSE)
file_identity <- function(path) {
  info <- file.info(path)
  list(path=normalizePath(path,mustWork=TRUE),sha256=hash(path),size=unname(info$size),
       mtime_utc=format(info$mtime,'%Y-%m-%dT%H:%M:%OS6Z',tz='UTC'))
}
paths <- list(image=manifest$copies$image$path,segmap=manifest$copies$segmap$path,
  original_image=manifest$originals_before$image$path,
  original_segmap=manifest$originals_before$segmap$path,
  v9='/home/uqpdas/WAVES_BD_pipeline/profuse_psf_v9.R',
  v13=file.path(root,'reference','psf_v13.R'),v14=file.path(root,'R','psf_v14.R'))
for(name in c('image','segmap')) {
  stopifnot(identical(hash(paths[[name]]),manifest$copies[[name]]$sha256),
    identical(verification$products[[name]]$numeric_sha256_f64le_R_columnmajor,
              manifest$copies[[name]]$array$numeric_sha256_f64le_R_columnmajor))
}
for(name in c('image','segmap')) stopifnot(
  identical(hash(paths[[paste0('original_',name)]]),manifest$originals_before[[name]]$sha256))
before <- lapply(paths,file_identity)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_json <- function(value,path) writeLines(jsonlite::toJSON(value,auto_unbox=TRUE,
  pretty=TRUE,digits=NA,null='null',na='null'),path)
dll <- getLoadedDLLs()[['ProFit']][['path']]
runtime <- list(R=R.version.string,job=Sys.getenv('SLURM_JOB_ID'),
  git_HEAD=system2('git',c('rev-parse','HEAD'),stdout=TRUE),package_manifest=package_manifest,
  ProFit_DLL=dll,ProFit_DLL_sha256=hash(dll),source_files_before=before,
  script_sha256=hash(file.path(root,'scripts','build_pilot_psf_candidate.R')),
  input_manifest_path=manifest_path,input_manifest_sha256=hash(manifest_path),
  input_Rfits_verification_sha256=hash(verification_path),
  scope='REQUESTED_PSF_STAR_FITS_ONLY_NO_OBSERVED_GALAXY_FIT_OR_SCIENTIFIC_GATE_DECISION')
write_json(runtime,file.path(out,'runtime.json'));print(runtime)
science <- bd_read_fits(paths$image,1L)
seg <- bd_read_fits(paths$segmap,1L)
grid <- bd_check_original_grid(science,seg)
position <- Rwcs::Rwcs_s2p(RA=manifest$ra_deg,Dec=manifest$dec_deg,keyvalues=seg$keyvalues)
stopifnot(is.matrix(position),all(c('x','y') %in% colnames(position)))
galaxy_xy <- as.numeric(position[1,c('x','y')])
anchor <- bd_segmap_converter(seg)(manifest$ra_deg,manifest$dec_deg)
stopifnot(seg$z[anchor$row,anchor$col]==4211L,
          identical(numeric_hash(science$z),manifest$copies$image$array$numeric_sha256_f64le_R_columnmajor),
          identical(numeric_hash(seg$z),manifest$copies$segmap$array$numeric_sha256_f64le_R_columnmajor))
image_identity <- object_hash(science$z)
seg_identity <- object_hash(seg$z)

matrix_measurement <- function(z) {
  if(!is.matrix(z) || any(!is.finite(z)) || any(z<0) || sum(z)<=0)
    return(list(status='INVALID_NUMERIC_KERNEL'))
  nr <- nrow(z);nc <- ncol(z)
  x <- seq_len(nr)-.5;y <- seq_len(nc)-.5
  centre <- c(sum(z*outer(x,rep(1,nc))),sum(z*outer(rep(1,nr),y)))/sum(z)
  radius <- sqrt(outer((x-nr/2)^2,(y-nc/2)^2,'+'))
  radii <- sort(unique(c(1,2,3,5,10,min(nr,nc)/4,min(nr,nc)/2)))
  list(status='NUMERIC_MEASUREMENT_ONLY',dimensions=dim(z),sum=sum(z),
    attributes=attributes(z),object_sha256=object_hash(z),numeric_sha256=numeric_hash(z),
    physical_centroid_xy=centre,physical_grid_centre_xy=c(nr,nc)/2,
    physical_centroid_offset_xy=centre-c(nr,nc)/2,
    normalized_180deg_asymmetry=max(abs(z-z[nr:1,nc:1,drop=FALSE]))/max(z),
    finite_kernel_encircled_fractions=data.frame(radius_pixels=radii,
      fraction=vapply(radii,function(r)sum(z[radius<=r])/sum(z),numeric(1))),
    finite_kernel_EE_interpretation='Fractions of already normalized truncated matrix; '
      %+% 'not fractions of infinite stellar flux and not measured stellar encircled energy.',
    builder_infinite_moffat_EE_approx=attr(z,'ee_infinite_moffat_approx',exact=TRUE),
    builder_infinite_moffat_lost_fraction_approx=attr(z,'ee_lost_fraction_approx',exact=TRUE),
    builder_EE_method=attr(z,'ee_method',exact=TRUE),
    builder_EE_interpretation='Analytic Moffat approximation, not independently measured stellar EE')
}
`%+%` <- function(a,b) paste0(a,b)
mutation_measurement <- function(psf) {
  baseline <- unserialize(serialize(psf,NULL))
  initial <- object_hash(psf)
  first <- NULL;maximum <- 0;first_change <- NA_integer_
  d <- pmax(as.integer(dim(psf))+20L,61L)
  for(i in seq_len(10L)) {
    invisible(ProFit::profitMakeModel(list(pointsource=list(xcen=d[1]/2,
      ycen=d[2]/2,mag=0)),psf=psf,magzero=0,dim=d,finesample=1L))
    maximum <- max(maximum,max(abs(psf-baseline)))
    current <- object_hash(psf)
    if(i==1L) first <- current
    if(is.na(first_change) && current!=initial) first_change <- i
  }
  list(scope='Ten synthetic point-source renders only; no F2 or observed galaxy fit',
    initial_sha256=initial,after_first_sha256=first,final_sha256=object_hash(psf),
    maximum_absolute_difference=maximum,first_changed_render=first_change,
    input_changed=!identical(initial,object_hash(psf)))
}

results <- list();summaries <- list();catalogue_hashes <- list();effective_settings <- list()
for(version in c('v13','v14')) {
  env <- new.env(parent=.GlobalEnv)
  sys.source(paths$v9,envir=env)
  sys.source(paths[[version]],envir=env)
  builder <- get('build_psf_for_image_v13',envir=env,inherits=FALSE)
  if(!identical(environment(builder),env)) stop('Builder environment isolation failed')
  # Deparse the pairlist as a whole, preserving required arguments as missing.
  defaults <- paste(deparse(formals(builder)),collapse='\n')
  settings <- list(image=science$z,segim=seg$z,magzero=0,
    galaxy_xy=galaxy_xy,exclude_segids=4211L,plot_profound=FALSE,plot_psf_fit=FALSE)
  if('build_v9_sidecars' %in% names(formals(builder))) settings$build_v9_sidecars <- FALSE
  effective_settings[[version]] <- list(builder_formal_defaults=defaults,
    explicit_arguments=list(image='EXACT_VALIDATED_IMAGE_ARRAY',segim='EXACT_ORIGINAL_LABELS',
      magzero=0,magzero_interpretation='COMPUTATIONAL_RELATIVE_UNITS_NOT_PHYSICAL_ZEROPOINT',
      galaxy_xy=galaxy_xy,exclude_segids=4211L,plot_profound=FALSE,plot_psf_fit=FALSE,
      build_v9_sidecars=settings$build_v9_sidecars),
    native_catalogue='ProFound static_photom=TRUE with original segim, sky=0, redosky=FALSE; no detection')
  set.seed(666L)
  started <- proc.time()[['elapsed']]
  warning_messages <- character();error_message <- NULL
  built <- tryCatch(withCallingHandlers(do.call(builder,settings),warning=function(w) {
    warning_messages <<- c(warning_messages,conditionMessage(w))
  }),error=function(e) {error_message <<- conditionMessage(e);NULL})
  elapsed <- proc.time()[['elapsed']]-started
  unchanged <- identical(object_hash(science$z),image_identity) && identical(object_hash(seg$z),seg_identity)
  saveRDS(built,file.path(out,paste0(version,'_builder_full.rds')),version=2,compress=TRUE)
  if(is.null(built)) {
    summaries[[version]] <- list(status='BUILDER_ERROR',error=error_message,
      warnings=warning_messages,elapsed_seconds=elapsed,source_inputs_unchanged=unchanged)
    next
  }
  results[[version]] <- built
  catalogue_hashes[[version]] <- object_hash(built$profound$segstats)
  candidates <- if(is.data.frame(built$candidates)) as.integer(built$candidates$segID) else integer()
  selected_ids <- if(length(built$selected_idx) && length(built$fits))
    vapply(built$fits[built$selected_idx],function(f)as.integer(f$segID),integer(1)) else integer()
  fits_summary <- lapply(built$fits,function(f) f[intersect(names(f),c('ok','segID','reason',
    'stamp_size','fit_par','fit_par_log','reduced_chisq','nllh','ndof','fit_limits',
    'opt_convergence','convergence','at_lim','sigma_pix','cutout_anchor',
    'segstats_xcen','segstats_ycen','segstats_mag','segstats_R50','segstats_R90','bg','bg_init'))])
  nondefault <- identical(built$output_kind,'gama') && is.matrix(built$psf_gama) &&
    length(selected_ids)>0L && is.numeric(attr(built$psf,'Nstar_used')) && attr(built$psf,'Nstar_used')>0
  numeric_kernel <- is.matrix(built$psf) && all(is.finite(built$psf)) &&
    all(built$psf>=0) && abs(sum(built$psf)-1)<1e-12
  summaries[[version]] <- list(status=if(!numeric_kernel) 'BLOCKED_INVALID_NUMERIC_KERNEL'
    else if(nondefault) 'RAW_MEASURED_PSF_CANDIDATE' else 'BLOCKED_DEFAULT_OR_ZERO_SELECTED_STARS',
    output_kind=built$output_kind,psf_quality_flag_raw=built$psf_quality_flag,
    warnings=warning_messages,elapsed_seconds=elapsed,source_inputs_unchanged=unchanged,
    star_catalogue_segstats_sha256=catalogue_hashes[[version]],candidate_star_ids=candidates,
    selected_star_ids=selected_ids,selected_indices=built$selected_idx,
    target4211_excluded_from_candidates=!4211L %in% candidates,
    star_fit_records=fits_summary,star_table=built$star_table,decision=built$decision,
    sky_pedestal_subtracted=built$sky_pedestal_subtracted,
    matrix_measurement=matrix_measurement(built$psf),
    mutation_measurement=if(numeric_kernel) mutation_measurement(built$psf) else NULL)
  if(is.data.frame(built$star_table)) write.table(built$star_table,
    file.path(out,paste0(version,'_star_table.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  if(is.data.frame(built$candidates)) write.table(built$candidates,
    file.path(out,paste0(version,'_candidate_stars.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  if(is.data.frame(built$profound$segstats)) saveRDS(built$profound$segstats,
    file.path(out,paste0(version,'_static_photometry_catalogue.rds')),version=2)
  write_json(summaries[[version]],file.path(out,paste0(version,'_measurements.json')))
}
source_after <- lapply(paths,file_identity)
sources_unchanged <- identical(before,source_after)
same_catalogue <- length(catalogue_hashes)==2L &&
  identical(catalogue_hashes$v13,catalogue_hashes$v14)
same_defaults <- identical(effective_settings$v13$builder_formal_defaults,
                           effective_settings$v14$builder_formal_defaults)
problems <- character()
if(!sources_unchanged) problems <- c(problems,'ORIGINAL_OR_SOURCE_FILE_CHANGED')
if(!same_catalogue) problems <- c(problems,'STATIC_STAR_CATALOGUES_DIFFER_OR_UNAVAILABLE')
if(!same_defaults) problems <- c(problems,'BUILDER_DEFAULT_SETTINGS_DIFFER')
for(version in names(summaries)) {
  s <- summaries[[version]]
  if(!isTRUE(s$source_inputs_unchanged)) problems <- c(problems,paste0(version,':RAW_ARRAY_MUTATION'))
  if(!identical(s$status,'RAW_MEASURED_PSF_CANDIDATE')) problems <- c(problems,paste0(version,':',s$status))
  if(!is.null(s$mutation_measurement) && isTRUE(s$mutation_measurement$input_changed))
    problems <- c(problems,paste0(version,':PSF_MUTATION'))
  if(!is.null(s$target4211_excluded_from_candidates) && !isTRUE(s$target4211_excluded_from_candidates))
    problems <- c(problems,paste0(version,':TARGET_IN_STAR_POOL'))
}
v14 <- summaries$v14
if(!is.null(results$v14) && identical(v14$status,'RAW_MEASURED_PSF_CANDIDATE') &&
   isTRUE(v14$source_inputs_unchanged) && !isTRUE(v14$mutation_measurement$input_changed) &&
   sources_unchanged && same_catalogue && same_defaults) {
  saveRDS(results$v14$psf,file.path(out,'v14_measured_psf_candidate.rds'),version=2)
}
report <- list(status=if(length(problems)) 'BLOCKED_RAW_PSF_DIAGNOSTIC' else 'RAW_PSF_CANDIDATE_AWAITING_ORIGINAL_CRITERIA',
  uberID='W11690092907692',runtime=runtime,source_files_after=source_after,
  source_files_unchanged=sources_unchanged,grid=grid,galaxy_original_WCS_xy=galaxy_xy,
  settings=effective_settings,static_star_catalogue_identical=same_catalogue,
  builder_defaults_identical=same_defaults,measurements=summaries,blockers=problems,
  physical_units_zeropoint_gain_certified=FALSE,magzero_interpretation='COMPUTATIONAL_RELATIVE_UNITS_ONLY',
  E1_E7_G5_status='ORIGINAL_CRITERIA_NOT_PROVIDED_NOT_EVALUATED',
  observed_galaxy_fits=0L,F2_per_object_runs=0L,config_modified=FALSE,approval_flags_modified=FALSE)
report$loaded_namespace_manifest <- lapply(sort(loadedNamespaces()),function(package) {
  desc <- utils::packageDescription(package)
  list(package=package,version=desc$Version,path=getNamespaceInfo(asNamespace(package),'path'),
       source_sha=desc$RemoteSha)
})
write_json(report,file.path(out,'psf_candidate_report.json'))
files <- list.files(out,full.names=TRUE)
writeLines(vapply(files,function(path)paste(hash(path),basename(path)),character(1)),
           file.path(out,'SHA256SUMS.txt'))
cat('RAW_PSF_CANDIDATE_MEASUREMENTS_RECORDED',report$status,'\n')
if(!sources_unchanged) stop('Source file mutation detected')
