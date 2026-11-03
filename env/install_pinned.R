#!/usr/bin/env Rscript
# NO historical Rlib writes. Always install into a new repo-local restart library.
source(file.path('R','00_load_all.R'))
root<-bd_find_root()
if(!requireNamespace('yaml',quietly=TRUE)) stop('Make yaml available to the active R; it reads the lock file')
if(!requireNamespace('remotes',quietly=TRUE)) stop('Make remotes available to the active R; no unpinned auto-bootstrap')
pins<-yaml::read_yaml(file.path(root,'config','package_pins.yaml'))$packages
lib<-file.path(root,'env','restart_Rlib')
if(file.exists(lib) && nzchar(Sys.readlink(lib))) stop('Refusing symlinked restart Rlib')
dir.create(lib,showWarnings=FALSE,recursive=TRUE)
if(any(normalizePath(.libPaths(),mustWork=FALSE)==normalizePath(lib,mustWork=TRUE)))
  message('Restart lib already active')
.libPaths(c(lib,.libPaths()))
# Install dependencies before consumers, without fetching unpinned Git remotes.
install_order<-c('ParmOff','celestial','magicaxis','LaplacesDemon','cmaeshpc',
                 'digest','Rfits','Rwcs','Highlander','ProFit','ProFound')
if(!setequal(install_order,names(pins))) stop('Update the pinned dependency order for this lock file')
for(n in install_order) {
  p<-pins[[n]]
  desc<-suppressWarnings(utils::packageDescription(n,lib.loc=lib))
  if(is.list(desc) && identical(desc$RemoteSha,p$sha)) {
    message(n,' already at requested SHA; skipping')
  } else {
    if(is.list(desc)) stop('Conflicting restart package exists: ',n, ' (will not overwrite it)')
    # remotes also searches other libraries; a matching historical install must
    # not prevent this missing package from being installed into the restart lib.
    remotes::install_github(paste0(p$repo,'@',p$sha),lib=lib,upgrade='never',
                            dependencies=FALSE,force=TRUE,build=TRUE,quiet=FALSE)
  }
  d<-suppressWarnings(utils::packageDescription(n,lib.loc=lib))
  if(!is.list(d) || !identical(d$RemoteSha,p$sha)) stop('Installed SHA mismatch for ',n)
  if(!is.null(p$version) && as.character(utils::packageVersion(n,lib.loc=lib))!=as.character(p$version))
    stop('Installed package version mismatch for ',n)
}
cat('EXACT_SOURCE_SHA_INSTALLS_VERIFIED\n')
