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
    saved_R_version=field(fit,'R.version'),
    psf_values=if(is.list(d)) d[['psf',exact=TRUE]] else NULL)
}
is_kernel <- function(z) is.numeric(z) && length(dim(z))==2L
# Physical-coordinate geometry of one kernel (pixel i centre = i-0.5; grid centre = dim/2).
kernel_geometry <- function(z) {
  if (!is_kernel(z) || !all(is.finite(z)) || !(sum(z)>0)) return(NULL)
  nr <- nrow(z);nc <- ncol(z);w <- z/sum(z)
  x <- seq_len(nr)-.5;y <- seq_len(nc)-.5
  cx <- sum(w*x);cy <- sum(w*rep(y,each=nr))
  list(dimensions=dim(z),sum=sum(z),min_value=min(z),peak_value=max(z),
    peak_index_xy=as.integer(arrayInd(which.max(z),dim(z))),
    physical_centroid_xy=c(cx,cy),centroid_offset_from_grid_centre_xy=c(cx-nr/2,cy-nc/2),
    second_moments=list(xx=sum(w*(x-cx)^2),yy=sum(w*rep((y-cy)^2,each=nr)),
      xy=sum(w*outer(x-cx,y-cy))),
    normalised_180deg_asymmetry=max(abs(z-z[nr:1,nc:1,drop=FALSE]))/max(abs(z)))
}
# Claude review P2-1: distinguish renormalisation/resizing from a different kernel.
# Same dimensions: EXACT (identical doubles), RENORMALISED (max|a/sum(a)-b/sum(b)| < tol),
# else DIFFERENT_KERNEL. Different dimensions: compare the centred common overlap,
# RESIZED when proportional there, else DIFFERENT_KERNEL_AND_DIMENSIONS.
classify_psf_values <- function(a,b,tol=1e-10) {
  if (!is_kernel(a) || !is_kernel(b)) return(list(class='NOT_COMPARABLE_NON_MATRIX'))
  if (!all(is.finite(a)) || !all(is.finite(b))) return(list(class='NOT_COMPARABLE_NON_FINITE'))
  sa <- sum(a);sb <- sum(b)
  if (!(sa>0) || !(sb>0)) return(list(class='NOT_COMPARABLE_NON_POSITIVE_SUM',sum_a=sa,sum_b=sb))
  da <- dim(a);db <- dim(b);common <- pmin(da,db)
  off_a <- (da-common)/2;off_b <- (db-common)/2
  if (any(off_a!=floor(off_a)) || any(off_b!=floor(off_b)))
    return(list(class='NOT_COMPARABLE_DIMENSION_PARITY',dims_a=da,dims_b=db))
  ca <- a[off_a[1]+seq_len(common[1]),off_a[2]+seq_len(common[2]),drop=FALSE]
  cb <- b[off_b[1]+seq_len(common[1]),off_b[2]+seq_len(common[2]),drop=FALSE]
  same_dims <- identical(as.integer(da),as.integer(db))
  na <- ca/sum(ca);nb <- cb/sum(cb)
  normalised_difference <- max(abs(na-nb))
  exact <- same_dims && identical(as.double(a),as.double(b))
  class <- if (exact) 'EXACT' else if (normalised_difference<tol) {
    if (same_dims) 'RENORMALISED' else 'RESIZED'
  } else if (same_dims) 'DIFFERENT_KERNEL' else 'DIFFERENT_KERNEL_AND_DIMENSIONS'
  list(class=class,tolerance=tol,same_dimensions=same_dims,dims_a=da,dims_b=db,
    common_dims=common,sum_a=sa,sum_b=sb,sum_ratio_a_over_b=sa/sb,
    max_abs_raw_difference_common=max(abs(ca-cb)),
    max_abs_normalised_difference=normalised_difference,
    normalised_difference_relative_to_peak=normalised_difference/max(abs(na),abs(nb)),
    pixels_differing_common=sum(ca!=cb),pixels_common=length(ca),
    flux_fraction_outside_common_a=1-sum(ca)/sa,flux_fraction_outside_common_b=1-sum(cb)/sb,
    centroid_offset_from_grid_centre_a=kernel_geometry(a)$centroid_offset_from_grid_centre_xy,
    centroid_offset_from_grid_centre_b=kernel_geometry(b)$centroid_offset_from_grid_centre_xy)
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
  accepted_psf <- accepted$psf_values;accepted$psf_values <- NULL
  source_psf <- source$psf_values;source$psf_values <- NULL
  cache_values <- list()
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
        if (is.numeric(p) && length(dim(p))==2L) {
          cache$kernels$ROOT <- psf_metadata(p);cache_values$ROOT <- p
        } else if (is.list(p)) for(key in names(p)) {
          v <- p[[key,exact=TRUE]]
          if(is.numeric(v)&&length(dim(v))==2L) {
            cache$kernels[[key]] <- psf_metadata(v)
            # Value comparison only for PSF kernels, not e.g. the cached segmentation map.
            if (grepl('^psf',key)) cache_values[[key]] <- v
          }
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
    accepted_psf_geometry=kernel_geometry(accepted_psf),
    source_psf_geometry=kernel_geometry(source_psf),
    accepted_source_psf_value_comparison=classify_psf_values(accepted_psf,source_psf),
    accepted_cache_kernel_value_comparisons=lapply(cache_values,
      function(p) classify_psf_values(accepted_psf,p)),
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
