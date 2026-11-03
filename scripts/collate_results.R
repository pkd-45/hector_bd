#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L || !args[[1]] %in% c('mac','friday')) stop('Usage: Rscript scripts/collate_results.R mac|friday RUN_ID')
source(file.path('R','00_load_all.R'));bd_load_all()
config<-bd_config(args[[1]]);print(bd_collate(config$paths$output_root,args[[2]]))
