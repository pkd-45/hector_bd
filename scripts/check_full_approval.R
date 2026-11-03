#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('profile needed')
source(file.path('R','00_load_all.R'));bd_load_all()
config<-bd_config(args[[1]])
if (!identical(config$pilot$full_832_approved,TRUE) ||
    !identical(config$pilot$pilot_approved,TRUE) ||
    !identical(config$pilot$pilot_review_completed,TRUE))
 stop('BLOCKED: full-run approval in config has not been granted')
bd_step0_gate(config)
bd_fingerprint(config)
cat('FULL_APPROVAL_FLAGS_PASS; source and pixel gates still enforced per galaxy\n')
