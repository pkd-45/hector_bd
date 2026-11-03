#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
# CLI is a thin adapter: same bd_run_one() and bd_stage() as IRkernel notebooks.
arg<-function(key,default=NULL) {
  hit<-grep(paste0('^--',key,'='),args,value=TRUE)
  if(length(hit)>1L) stop('Duplicate option ',key)
  if(!length(hit)) return(default)
  substring(hit, nchar(key)+4L)
}
scriptarg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(length(scriptarg)!=1L) stop('Unable to identify CLI script root')
root <- normalizePath(file.path(dirname(sub('^--file=','',scriptarg)), '..'),mustWork=TRUE)
if(!file.exists(file.path(root,'.hector_bd_root'))) stop('Repository marker absent')
setwd(root)
source(file.path(root,'R','00_load_all.R'))
.libPaths(c(file.path(root,'env','restart_Rlib'),.libPaths()))
bd_load_all(root)
profile<-arg('profile','mac');cfg<-bd_config(profile,root)
id<-arg('id')
if(is.null(id)) stop('Provide --id=W...')
mode<-arg('mode','synthetic');run_id<-arg('run-id','test_synthetic')
stages<-if(mode=='synthetic') c('inputs','g2','f1_f9') else c('inputs','g2','f1_f9','fit','postprocess','card')
if(!is.null(arg('stages'))) {
  s<-strsplit(arg('stages'),',',fixed=TRUE)[[1]]
  if(mode!='synthetic' || !identical(s,c('inputs','g2'))) stop('Stage override permitted ONLY for synthetic inputs,g2 test')
  stages<-s
}
bd_run_one(id,cfg,run_id,mode,
  output_override=if(mode=='synthetic' && nzchar(Sys.getenv('BD_OUTPUT_OVERRIDE'))) Sys.getenv('BD_OUTPUT_OVERRIDE') else NULL,
  stages=stages)
cat('Shared runner completed stage sequence:',paste(stages,collapse=' -> '),'\n')
