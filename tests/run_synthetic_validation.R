#!/usr/bin/env Rscript
# Execute the existing recovery assertions and preserve diagnostics on failure.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('Synthetic validation requires SLURM')
out <- Sys.getenv('BD_SYNTHETIC_REPORT_DIR')
if (!nzchar(out) || !dir.exists(out)) stop('Missing synthetic report directory')
source('R/00_load_all.R');bd_load_all()
stopifnot(requireNamespace('ProFit',quietly=TRUE),requireNamespace('digest',quietly=TRUE),
  requireNamespace('jsonlite',quietly=TRUE))
lock <- yaml::read_yaml('config/package_pins.yaml')$packages$ProFit
desc <- utils::packageDescription('ProFit')
stopifnot(identical(desc$RemoteSha,lock$sha),
  as.character(utils::packageVersion('ProFit'))==as.character(lock$version))
ns_path <- normalizePath(getNamespaceInfo(asNamespace('ProFit'),'path'))
expected <- normalizePath(file.path(Sys.getenv('R_LIBS_USER'),'ProFit'))
stopifnot(identical(ns_path,expected))
dll <- getLoadedDLLs()[['ProFit']][['path']]
runtime <- list(R=R.version.string,ProFit_version=desc$Version,
  ProFit_sha=desc$RemoteSha,ProFit_path=ns_path,ProFit_DLL=dll,
  ProFit_DLL_sha256=digest::digest(file=dll,algo='sha256'),
  git_HEAD=system2('git',c('rev-parse','HEAD'),stdout=TRUE),
  source_sha256=bd_source_hash(),R_LIBS_USER=Sys.getenv('R_LIBS_USER'),
  job=Sys.getenv('SLURM_JOB_ID'))
writeLines(jsonlite::toJSON(runtime,auto_unbox=TRUE,pretty=TRUE),
  file.path(out,'runtime.json'))
print(runtime)
test <- new.env(parent=.GlobalEnv)
failure <- tryCatch({source('tests/test_synthetic_recovery.R',local=test);NULL},
  error=function(e) conditionMessage(e))
report <- list(status=if(is.null(failure)) 'PASS' else 'FAIL',error=failure,
  runtime=runtime,elapsed_seconds=if(exists('elapsed',test,inherits=FALSE)) test$elapsed else NULL)
if (exists('truth',test,inherits=FALSE)) report$truth <- test$truth
if (exists('fits',test,inherits=FALSE)) {
  report$fits <- lapply(test$fits,function(f) f[c('objective','parameters',
    'modellist_best','chosen_start','finesample','psf_type','psf_renderer_flag',
    'convergence','convergence_message','evaluations_total','failed_evaluations',
    'uncertified_starts','starts','psf_unchanged','status')])
  fixed <- test$fits$D_FIXED4$modellist_best$sersic
  report$physical_recovery <- data.frame(
    parameter=c('xcen','ycen','bulge_mag','disc_mag','bulge_re','disc_re','disc_axrat','disc_ang'),
    truth=c(28.7,33.3,18.6,17.9,3.5,14,.55,32),
    recovered=c(fixed$xcen[1],fixed$ycen[1],fixed$mag[1],fixed$mag[2],
      fixed$re[1],fixed$re[2],fixed$axrat[2],fixed$ang[2]),
    tolerance=c(1,1,.65,.65,1.5,4,.18,15))
  r <- report$physical_recovery
  r$absolute_error <- abs(r$recovered-r$truth)
  r$absolute_error[8] <- abs(((r$recovered[8]-r$truth[8]+90) %% 180)-90)
  r$pass <- r$absolute_error<r$tolerance
  report$physical_recovery <- r
  print(r)
}
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,null='null',
  na='null',digits=NA),file.path(out,'synthetic_recovery_report.json'))
if (!is.null(failure)) stop('SYNTHETIC_RECOVERY_FAILED: ',failure)
cat('SYNTHETIC_VALIDATION_REPORT_PASS\n')
