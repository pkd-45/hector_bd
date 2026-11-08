#!/usr/bin/env Rscript
# Verify the installed lock and exercise native code only inside a SLURM job.
if (!nzchar(Sys.getenv('SLURM_JOB_ID'))) stop('Run verification through SLURM')
source(file.path('R','00_load_all.R'))
root <- bd_find_root()
lib <- normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
lock <- yaml::read_yaml(file.path(root,'config','package_pins.yaml'))
report <- Sys.getenv('BD_ENV_REPORT_DIR',unset=file.path(root,'outputs',
  paste0('r_verification_',Sys.getenv('SLURM_JOB_ID'))))
dir.create(report,recursive=TRUE,showWarnings=FALSE)
cat('R:',R.version.string,'\nLibrary:',lib,'\n')
packages <- c(names(lock$packages),names(lock$cran_runtime))
manifest <- lapply(packages,function(package) {
  path <- normalizePath(file.path(lib,package),mustWork=TRUE)
  if (!identical(dirname(path),lib)) stop('Package outside restart library: ',package)
  d <- utils::packageDescription(package,lib.loc=lib)
  if (!is.list(d)) stop('Missing package: ',package)
  p <- lock$packages[[package]]
  expected_version <- if (is.null(p)) lock$cran_runtime[[package]] else p$version
  if (!is.null(expected_version) &&
      utils::packageVersion(package,lib.loc=lib) != package_version(as.character(expected_version)))
    stop('Version mismatch: ',package)
  if (!is.null(p) && (!identical(d$RemoteSha,p$sha) ||
      !identical(paste(d$RemoteUsername,d$RemoteRepo,sep='/'),p$repo)))
    stop('Git source mismatch: ',package)
  ns <- loadNamespace(package,lib.loc=lib)
  if (!identical(normalizePath(getNamespaceInfo(ns,'path')),path))
    stop('Loaded namespace outside restart library: ',package)
  data.frame(package=package,version=as.character(utils::packageVersion(package,lib.loc=lib)),
    source=if(is.null(p)) 'CRAN' else p$repo,
    source_sha=if(is.null(p)) '' else d$RemoteSha,path=path,
    needs_compilation=if(is.null(d$NeedsCompilation)) '' else d$NeedsCompilation,
    stringsAsFactors=FALSE)
})
manifest <- do.call(rbind,manifest)
write.table(manifest,file.path(report,'package_manifest.tsv'),sep='\t',
  row.names=FALSE,quote=FALSE)
print(manifest[,c('package','version','source_sha')],row.names=FALSE)

# Exercise both bundled cubature archives, the maps nested build, and Fortran.
h <- cubature::hcubature(function(x) sum(x*x),c(0,0),c(1,1),tol=1e-8)
c <- cubature::cuhre(function(x) sum(x*x),lowerLimit=c(0,0),
  upperLimit=c(1,1),relTol=1e-8)
stopifnot(abs(h$integral-2/3)<1e-7,abs(c$integral-2/3)<1e-7)
world <- maps::map('world',plot=FALSE)
stopifnot(length(world$x)>100L,any(is.finite(world$x)))
projected <- mapproj::mapproject(c(0,1),c(0,1),projection='mercator')
stopifnot(all(is.finite(projected$x)),all(is.finite(projected$y)))
nearest <- RANN::nn2(cbind(0:3,0:3),matrix(c(0,0),nrow=1),k=1L)
stopifnot(nearest$nn.idx[1,1]==1L,nearest$nn.dists[1,1]==0)
stopifnot(identical(matrixStats::rowSums2(matrix(1:6,nrow=2)),c(9,12)))
density <- sm::sm.density(seq(-2,2,length.out=20),display='none')
stopifnot(all(is.finite(density$estimate)),all(density$estimate>=0))
dt <- data.table::data.table(x=1:4,g=c(1L,1L,2L,2L))
stopifnot(identical(dt[,sum(x),by=g]$V1,c(3L,7L)))
stopifnot(identical(digest::digest('abc',algo='sha256',serialize=FALSE),
  'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'))
cat('COMPILED_RUNTIME_SMOKE_TESTS_PASS\n')

dlls <- getLoadedDLLs()
dll_paths <- unique(vapply(dlls,function(d) d[['path']],character(1)))
dll_paths <- dll_paths[file.exists(dll_paths)]
# Running R already loads libR globally. A standalone ldd does not, and Conda
# deliberately avoids exporting LD_LIBRARY_PATH from R's Linux startup script.
# Give only the diagnostic subprocess R's own library directories; namespace
# loading above verifies the real runtime without altering its environment.
link_search <- paste(Filter(nzchar,c(R.home('lib'),dirname(R.home()),
                                    Sys.getenv('LD_LIBRARY_PATH'))),collapse=':')
cat('LDD_SEARCH_PATH',link_search,'\n')
native_manifest <- lapply(dll_paths,function(path) {
  path <- normalizePath(path,mustWork=TRUE)
  links <- system2('ldd',shQuote(path),stdout=TRUE,stderr=TRUE,
                   env=paste0('LD_LIBRARY_PATH=',shQuote(link_search)))
  cat('\nDLL:',path,'\n',paste(links,collapse='\n'),'\n',
    file=file.path(report,'shared_library_links.txt'),append=TRUE)
  if ((!is.null(attr(links,'status')) && attr(links,'status')!=0L) ||
      any(grepl('not found',links,fixed=TRUE)))
    stop('Unresolved DLL dependency: ',path,'\n',paste(links,collapse='\n'))
  data.frame(path=path,sha256=digest::digest(file=path,algo='sha256'),
    stringsAsFactors=FALSE)
})
native_manifest <- do.call(rbind,native_manifest)
write.table(native_manifest,file.path(report,'shared_library_manifest.tsv'),
  sep='\t',row.names=FALSE,quote=FALSE)
for (package in packages) {
  native <- list.files(file.path(lib,package,'libs'),pattern='\\.so$',
    recursive=TRUE,full.names=TRUE)
  if (length(native) && !all(normalizePath(native) %in% native_manifest$path))
    stop('Installed shared library was not loaded: ',package)
}
namespaces <- lapply(loadedNamespaces(),function(package) {
  ns <- asNamespace(package)
  path <- if (identical(package,'base')) file.path(.Library,'base') else
    getNamespaceInfo(ns,'path')
  d <- utils::packageDescription(package,lib.loc=dirname(path))
  data.frame(package=package,version=as.character(utils::packageVersion(package,lib.loc=dirname(path))),
    source_sha=if(is.null(d$RemoteSha)) '' else d$RemoteSha,
    path=path,stringsAsFactors=FALSE)
})
write.table(do.call(rbind,namespaces),file.path(report,'loaded_namespace_manifest.tsv'),
  sep='\t',row.names=FALSE,quote=FALSE)
cat('PINNED_NAMESPACES_AND_SHARED_LIBRARIES_VERIFIED\n')
