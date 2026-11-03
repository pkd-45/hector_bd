#!/usr/bin/env Rscript
# D1 read-only extractor. No package imports, no model execution, no fitting.
# VERIFIED_FROM_SOURCE_REVIEW: historical inputs.rds and fit_*.rds are R serialized.
# REPORTED: fit data may live at fit$Data, Data or nested fit$fit$Data.
# INFERRED/UNRESOLVED: object fields differ between versions; discover and log.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop('usage: export_rds_arrays.R <inputs.rds> <working-output-folder> <fit-files-separated-by-|>')
in_path <- args[[1]]; out <- args[[2]]; fit_args <- strsplit(args[[3]], '|', fixed=TRUE)[[1]]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
if (!dir.exists(out)) stop('cannot create output')
metadata <- list(); arrays <- list(); notes <- character()
sanitize <- function(x) {
  x <- paste(x, collapse=',')
  x <- gsub('[\t\r\n]', ' ', x)
  substr(x, 1L, 4000L)
}
append_meta <- function(src, key, val, evidence='OBSERVED_RDS_FIELD') {
  if (is.null(val)) return(invisible(NULL))
  if (is.factor(val)) val <- as.character(val)
  if (length(val)>20L || (is.list(val) && length(val)>3L)) return(invisible(NULL))
  if (is.atomic(val) && is.null(dim(val))) {
    metadata[[length(metadata)+1L]] <<- data.frame(source=src, key=key, value=sanitize(val), class=paste(class(val),collapse=','), evidence=evidence, stringsAsFactors=FALSE)
  }
}
export_matrix <- function(src, key, obj) {
  if (is.null(obj) || length(dim(obj)) != 2L || !(is.numeric(obj) || is.logical(obj))) return(FALSE)
  if (!all(dim(obj) > 0L)) return(FALSE)
  filename <- paste0(gsub('[^[:alnum:]_]+','_',paste(src,key,sep='__')), '.f64')
  for (at in c('fwhm','fwhm_arcsec','beta','con','axrat','ang','Nstar_used','size_px','target_flux_frac','encircled_energy')) {
    av <- attr(obj,at,exact=TRUE)
    if (!is.null(av)) append_meta(src,paste0('attr.',key,'.',at),av,'VERIFIED_RDS_ATTRIBUTE')
  }
  con <- file(file.path(out,filename), 'wb')
  on.exit(close(con),add=TRUE)
  writeBin(as.double(obj), con, size=8L, endian='little')
  arrays[[length(arrays)+1L]] <<- data.frame(source=src, key=key, file=filename, nrow=nrow(obj), ncol=ncol(obj), class=paste(class(obj),collapse=','), stringsAsFactors=FALSE)
  TRUE
}
visit_scalars <- function(x,src,root='',depth=0L) {
  if (!is.list(x) || depth>3L) return(invisible(NULL))
  keys <- names(x)
  if (is.null(keys)) return(invisible(NULL))
  for (key in keys) {
    if (!nzchar(key)) next
    v <- tryCatch(x[[key]],error=function(e) NULL)
    child <- if (nzchar(root)) paste(root,key,sep='.') else key
    if (is.atomic(v) && is.null(dim(v))) append_meta(src,child,v)
    else if (is.list(v) && length(v)<=12L && !is.data.frame(v)) visit_scalars(v,src,child,depth+1L)
  }
}
read_safe <- function(path,src) {
  if (!file.exists(path)) {notes <<- c(notes,paste(src,'MISSING_FILE',path,sep='\t'));return(NULL)}
  tryCatch(readRDS(path), error=function(e) {notes <<- c(notes,paste(src,'READ_ERROR',sanitize(conditionMessage(e)),sep='\t'));NULL})
}
locate_data <- function(obj,depth=0L) {
  if (depth>4L || !is.list(obj)) return(NULL)
  if (is.list(obj$Data) && (!is.null(obj$Data$region) || !is.null(obj$Data$image))) return(obj$Data)
  if (!is.null(obj$region) && !is.null(obj$image)) return(obj)
  for (key in c('fit','output','result','results','res','best','accepted','object')) {
    if (!is.null(obj[[key]])) {
      a <- locate_data(obj[[key]],depth+1L)
      if (!is.null(a)) return(a)
    }
  }
  NULL
}
locate_profound <- function(obj, depth=0L) {
  if (depth>4L || !is.list(obj)) return(NULL)
  if (is.list(obj$profound) && any(c('image','segim','segim_orig') %in% names(obj$profound))) return(obj$profound)
  for (key in c('fit','output','result','results','res','best','accepted','object')) {
    if (!is.null(obj[[key]])) {
      z <- locate_profound(obj[[key]], depth+1L)
      if (!is.null(z)) return(z)
    }
  }
  NULL
}
inputs <- read_safe(in_path,'inputs')
if (is.list(inputs)) {
  append_meta('inputs','__names__',paste(names(inputs),collapse=','))
  visit_scalars(inputs,'inputs')
  for (key in c('image','segim','sigma','mask','region')) export_matrix('inputs',key,inputs[[key]])
} else notes <- c(notes,'inputs\tNOT_A_LIST_OR_UNREADABLE')
for (i in seq_along(fit_args)) {
  path <- fit_args[[i]]
  if (!nzchar(path)) next
  nm <- tools::file_path_sans_ext(basename(path))
  o <- read_safe(path,nm)
  if (is.null(o)) next
  if (is.list(o)) {append_meta(nm,'__names__',paste(names(o),collapse=',')); visit_scalars(o,nm)}
  profound <- locate_profound(o)
  for (key in c('image','segim','segim_orig','mask')) {
    saved <- !is.null(profound) && !is.null(profound[[key]])
    append_meta(nm,paste0('profound_presence_',key),if (saved) 'PRESENT' else 'ABSENT','VERIFIED_RDS_FIELD_PRESENCE')
    if (saved) export_matrix(nm,paste0('profound_',key),profound[[key]])
  }
  if (!is.null(profound)) append_meta(nm,'__profound_names__',paste(names(profound),collapse=','))
  d <- locate_data(o)
  if (is.null(d)) {notes <- c(notes,paste(nm,'DATA_NOT_FOUND',path,sep='\t'));next}
  append_meta(nm,'__Data_names__',paste(names(d),collapse=','))
  visit_scalars(d,nm,'Data')
  # Historical rendering-branch evidence; absent fields are explicitly recorded.
  append_meta(nm,'Data_finesample',if(is.null(d$finesample)) 'ABSENT' else d$finesample,'VERIFIED_RDS_FIELD_PRESENCE')
  append_meta(nm,'Data_psftype',if(is.null(d$psftype)) 'ABSENT' else d$psftype,'VERIFIED_RDS_FIELD_PRESENCE')
  for (key in c('image','segim','sigma','mask','region','calcregion','psf')) export_matrix(nm,key,d[[key]])
}
psfp <- file.path(dirname(in_path),'psf.rds')
p <- read_safe(psfp,'psf_rds')
if (is.list(p)) {
  append_meta('psf_rds','__names__',paste(names(p),collapse=','))
  visit_scalars(p,'psf_rds')
  for (key in c('psf','psf_gama','psf_default','psf_empirical','psf_apodized')) export_matrix('psf_rds',key,p[[key]])
}
empty_meta <- data.frame(source=character(),key=character(),value=character(),class=character(),evidence=character())
empty_arrays <- data.frame(source=character(),key=character(),file=character(),nrow=integer(),ncol=integer(),class=character())
write.table(if (length(metadata)) do.call(rbind,metadata) else empty_meta,file.path(out,'metadata.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='')
write.table(if (length(arrays)) do.call(rbind,arrays) else empty_arrays,file.path(out,'arrays.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='')
writeLines(notes,file.path(out,'notes.txt'))
cat('RDS_EXPORT=COMPLETE\n')
