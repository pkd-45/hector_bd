#!/usr/bin/env Rscript
# Read ONLY newly decompressed copies with the actual pinned pipeline adapter.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('Actual FITS reads require SLURM')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L || file.exists(args[2])) stop('Expected manifest and NEW verification output')
source('R/00_load_all.R')
root <- bd_find_root()
lib <- normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace('jsonlite',quietly=TRUE),requireNamespace('yaml',quietly=TRUE),
          requireNamespace('digest',quietly=TRUE))
pins <- yaml::read_yaml(file.path(root,'config','package_pins.yaml'))$packages
for(package in c('Rfits','Rwcs','digest')) {
  desc <- utils::packageDescription(package,lib.loc=lib)
  stopifnot(identical(desc$RemoteSha,pins[[package]]$sha))
  loadNamespace(package,lib.loc=lib)
  stopifnot(identical(normalizePath(getNamespaceInfo(asNamespace(package),'path')),
                    normalizePath(file.path(lib,package))))
}
source('R/config_io.R');source('R/io_fits.R')
manifest <- jsonlite::fromJSON(args[1],simplifyVector=FALSE)
stopifnot(identical(manifest$uberID,'W11690092907692'))
products <- list();read <- list()
for(name in c('image','segmap')) {
  path <- manifest$copies[[name]]$path
  stopifnot(!grepl('[.]gz$',path),startsWith(normalizePath(path),paste0(root,'/outputs/')))
  obj <- bd_read_fits(path,1L)
  z <- obj$z
  evidence <- manifest$copies[[name]]$array
  samples <- lapply(evidence$samples,function(p) list(
    fits_x_1based=p$fits_x_1based,fits_y_1based=p$fits_y_1based,
    value=unname(z[p$fits_x_1based,p$fits_y_1based]),
    finite=is.finite(z[p$fits_x_1based,p$fits_y_1based])))
  products[[name]] <- list(matrix_dims=dim(z),finite_count=sum(is.finite(z)),
    copied_fits_sha256=bd_file_hash(path),samples=samples,
    numeric_sha256_f64le_R_columnmajor=digest::digest(
      writeBin(as.double(z),raw(),size=8L,endian='little'),algo='sha256',serialize=FALSE))
  read[[name]] <- obj
}
grid <- bd_check_original_grid(read$image,read$segmap)
position <- bd_segmap_converter(read$segmap)(manifest$ra_deg,manifest$dec_deg)
label <- read$segmap$z[position$row,position$col]
stopifnot(as.integer(label)==4211L)
runtime <- list(R=R.version.string,R_LIBS_USER=Sys.getenv('R_LIBS_USER'),
  Rfits_version=as.character(utils::packageVersion('Rfits')),
  Rfits_sha=utils::packageDescription('Rfits')$RemoteSha,
  Rwcs_version=as.character(utils::packageVersion('Rwcs')),
  Rwcs_sha=utils::packageDescription('Rwcs')$RemoteSha,
  scope='COPY_READ_ORIENTATION_ONLY_NO_SCIENTIFIC_GATE_APPROVAL')
result <- list(runtime=runtime,products=products,grid=grid,
  target_R_1based_position=position,target_label=label)
writeLines(jsonlite::toJSON(result,auto_unbox=TRUE,pretty=TRUE,digits=NA,
  null='null',na='null'),args[2])
cat('RFITS_EXPLICIT_HDU_AND_ORIENTATION_EVIDENCE_RECORDED\n')
