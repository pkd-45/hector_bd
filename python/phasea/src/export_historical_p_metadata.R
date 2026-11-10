#!/usr/bin/env Rscript
# Read historical serialized fields only: no restore helpers, models or fits.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('Bulk RDS export requires SLURM')
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop('usage: export_historical_p_metadata.R bindings.tsv output_dir')
if (!requireNamespace('digest',quietly=TRUE) || !requireNamespace('jsonlite',quietly=TRUE))
  stop('digest and jsonlite must be available')
binding <- read.delim(args[[1]],stringsAsFactors=FALSE,colClasses='character',
  check.names=FALSE,quote='"',comment.char='')
if (nrow(binding)!=68L || anyDuplicated(binding$uberID)) stop('Expected exact68 unique P bindings')
out <- normalizePath(args[[2]],mustWork=TRUE)
destination <- file.path(out,'historical_p_metadata.jsonl')
if (file.exists(destination)) stop('Refusing metadata overwrite')
json <- function(x) jsonlite::toJSON(x,auto_unbox=TRUE,null='null',na='null',digits=NA)
hash_file <- function(path) digest::digest(file=path,algo='sha256')
small <- function(x,depth=0L) {
  if (is.null(x)) return(NULL)
  if (inherits(x,c('numeric_version','package_version'))) return(as.character(x))
  if (is.atomic(x) && is.null(dim(x)) && length(x)<=100L) return(unname(x))
  if (is.list(x) && depth<4L && length(x)<=30L) return(lapply(x,small,depth=depth+1L))
  list(class=class(x),length=length(x),dimensions=dim(x),value='NOT_EXPANDED')
}
field <- function(x,key) {
  if (!is.list(x) || !key %in% names(x)) return(list(presence='ABSENT',value=NULL))
  v <- x[[key,exact=TRUE]]
  list(presence=if(is.null(v)) 'PRESENT_NULL' else 'PRESENT',value=small(v))
}
psf_metadata <- function(p) {
  if (is.null(p)) return(list(presence='ABSENT',value=NULL))
  ans <- list(presence='PRESENT',class=class(p),dimensions=dim(p),
    attributes=small(attributes(p)))
  if (is.numeric(p) && length(dim(p))==2L) {
    ans$numeric_sha256_f64le_columnmajor <- digest::digest(
      writeBin(as.double(p),raw(),size=8L,endian='little'),algo='sha256',serialize=FALSE)
    ans$serialized_sha256_v2 <- digest::digest(serialize(p,NULL,version=2),
      algo='sha256',serialize=FALSE)
    ans$sum <- sum(p)
    ans$finite_count <- sum(is.finite(p))
    ans$pixel_count <- length(p)
  } else ans$numeric_sha256_f64le_columnmajor <- NULL
  ans
}
extract_fit <- function(path) {
  obj <- readRDS(path)
  fit <- if(is.list(obj)) obj[['fit',exact=TRUE]] else NULL
  d <- if(is.list(fit)) fit[['Data',exact=TRUE]] else NULL
  list(path=path,file_sha256=hash_file(path),top_level_names=names(obj),
    fit_names=names(fit),Data_names=names(d),
    Data_finesample=field(d,'finesample'),Data_psftype=field(d,'psftype'),
    Data_psffinesampled=field(d,'psffinesampled'),
    Data_magzero=field(d,'magzero'),Data_like_func=field(d,'like.func'),
    Data_psf=psf_metadata(if(is.list(d)) d[['psf',exact=TRUE]] else NULL),
    Data_modellist_psf=field(if(is.list(d)) d[['modellist',exact=TRUE]] else NULL,'psf'),
    final_modellist_psf=field(if(is.list(fit)) fit$finalmodel$modellist else NULL,'psf'),
    saved_ProFit_version=field(fit,'ProFit.version'),
    saved_ProFound_version=field(fit,'ProFound.version'),
    saved_Highlander_version=field(fit,'Highlander.version'),
    saved_ProFit_sha=field(fit,'ProFit.sha'),
    saved_R_version=field(fit,'R.version'))
}
compare_psf <- function(a,b) {
  if (!identical(a$presence,'PRESENT') || !identical(b$presence,'PRESENT')) return('MISSING_PSF_METADATA')
  if (is.null(a$numeric_sha256_f64le_columnmajor) || is.null(b$numeric_sha256_f64le_columnmajor))
    return('NON_MATRIX_PSF_NOT_COMPARED')
  if (!identical(a$dimensions,b$dimensions)) return('DIFFERENT_DIMENSIONS')
  if (identical(a$numeric_sha256_f64le_columnmajor,b$numeric_sha256_f64le_columnmajor))
    'EXACT_NUMERIC_KERNEL_MATCH' else 'DIFFERENT_NUMERIC_KERNEL_VALUES'
}
runtime <- list(R=R.version.string,job=Sys.getenv('SLURM_JOB_ID'),
  digest_version=as.character(utils::packageVersion('digest')),
  digest_path=find.package('digest'),scope='EXTRACTION_RUNTIME_NOT_HISTORICAL_RUNTIME_PROOF')
writeLines(json(runtime),file.path(out,'metadata_extraction_runtime.json'))
con <- file(destination,'wt')
for (i in seq_len(nrow(binding))) {
  b <- as.list(binding[i,,drop=FALSE])
  accepted <- extract_fit(b$accepted_fit_rds)
  source <- extract_fit(b$source_rds)
  flags <- character()
  for (label in c('accepted','source')) {
    fit <- if(label=='accepted') accepted else source
    for (name in c('Data_finesample','Data_psftype','Data_psffinesampled',
                  'saved_ProFit_version','saved_ProFit_sha','saved_R_version'))
      if (!identical(fit[[name]]$presence,'PRESENT')) flags <- c(flags,paste(label,name,fit[[name]]$presence,sep=':'))
    if (!identical(fit$Data_psf$presence,'PRESENT')) flags <- c(flags,paste0(label,':Data_psf:ABSENT'))
  }
  cache_path <- b$psf_recorded_cache_path
  cache <- list(path=cache_path,path_evidence=b$psf_cache_path_status,
    status='NO_UNIQUE_EXPLICIT_LOG_REFERENCE',kernels=list())
  if (nzchar(cache_path)) {
    if (file.exists(cache_path)) {
      cache_read <- tryCatch(list(ok=TRUE,value=readRDS(cache_path)),
        error=function(e) list(ok=FALSE,error=conditionMessage(e)))
      if (!cache_read$ok) {
        cache$status <- 'EXPLICIT_LOG_REFERENCE_READ_ERROR'
        cache$read_error <- cache_read$error
      } else {
        p <- cache_read$value
        cache$status <- 'EXPLICIT_LOG_REFERENCE_READ'
        cache$file_sha256 <- hash_file(cache_path)
        cache$names <- names(p)
        cache$scalar_metadata <- if(is.list(p)) lapply(p,function(v)
          if(is.atomic(v)&&is.null(dim(v))&&length(v)<=100L) small(v) else NULL) else NULL
        if (is.numeric(p) && length(dim(p))==2L) cache$kernels$ROOT <- psf_metadata(p)
        else if (is.list(p)) for(key in names(p)) {
          v <- p[[key,exact=TRUE]]
          if(is.numeric(v)&&length(dim(v))==2L) cache$kernels[[key]] <- psf_metadata(v)
        }
      }
    } else cache$status <- 'EXPLICIT_LOG_REFERENCE_FILE_MISSING'
  }
  if (!identical(cache$status,'EXPLICIT_LOG_REFERENCE_READ')) flags <- c(flags,paste0('PSF_CACHE:',cache$status))
  comparisons <- lapply(cache$kernels,function(p) compare_psf(accepted$Data_psf,p))
  record <- list(uberID=b$uberID,model='P',selection_layer='R61',fit_layer='R53',
    run_index=b$run_index,start_id=b$start_id,
    accepted=accepted,source=source,
    accepted_source_psf_comparison=compare_psf(accepted$Data_psf,source$Data_psf),
    psf_cache=cache,accepted_cache_kernel_comparisons=comparisons,
    psf_log_evidence=jsonlite::fromJSON(b$psf_log_lines_json,simplifyVector=FALSE),
    source_run_log=b$source_run_log,source_run_log_sha256=b$source_run_log_sha256,
    metadata_flags=I(unique(flags)),
    interpretation='ABSENT_FIELDS_ARE_UNKNOWN_NO_HISTORICAL_SETTING_OR_DEFECT_INFERRED')
  writeLines(json(record),con)
  flush(con)
  cat(sprintf('[%d/68] %s run=%s metadata_flags=%d\n',i,b$uberID,b$run_index,length(flags)))
}
close(con)
cat('HISTORICAL_P_RDS_METADATA_EXPORT_COMPLETE\n')
