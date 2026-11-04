#!/usr/bin/env Rscript
# Uses the active R with the separate pinned restart library.
source(file.path('R','00_load_all.R'))
root<-bd_find_root()
lib<-normalizePath(file.path(root,'env','restart_Rlib'),mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
if(!requireNamespace('IRkernel',quietly=TRUE))
  stop('Install a loadable IRkernel and its dependencies for the active R; see the Mac setup in README.md')

name<-'waves_r_bd'
IRkernel::installspec(name=name,displayname='R (waves_r BD pinned)',user=TRUE)

# Jupyter starts kernels independently of the shell that registered them. Persist
# the restart library in the kernelspec instead of relying on inherited R_LIBS_USER.
kernel_list<-system2('jupyter',c('kernelspec','list','--json'),stdout=TRUE,stderr=TRUE)
status<-attr(kernel_list,'status')
if(!is.null(status) && status!=0) stop('Unable to locate the registered Jupyter kernelspec')
kernels<-jsonlite::fromJSON(paste(kernel_list,collapse='\n'),simplifyVector=FALSE)
resource_dir<-kernels$kernels[[name]]$resource_dir
if(is.null(resource_dir) || !nzchar(resource_dir)) stop('Registered kernelspec not found: ',name)
spec_file<-file.path(resource_dir,'kernel.json')
spec<-jsonlite::fromJSON(spec_file,simplifyVector=FALSE)
existing<-strsplit(Sys.getenv('R_LIBS_USER'),.Platform$path.sep,fixed=TRUE)[[1]]
spec$env$R_LIBS_USER<-paste(unique(c(lib,existing[nzchar(existing)])),collapse=.Platform$path.sep)
spec$argv<-c(spec$argv[1:2],
             '-e',
             paste0('.libPaths(c(',dQuote(lib),',.libPaths()));IRkernel::main()'),
             spec$argv[(match('-e',spec$argv)+2):length(spec$argv)])
jsonlite::write_json(spec,spec_file,auto_unbox=TRUE,pretty=TRUE)
cat('Registered ',name,' with R_LIBS_USER=',spec$env$R_LIBS_USER,'\n',sep='')
