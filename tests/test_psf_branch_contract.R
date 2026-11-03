#!/usr/bin/env Rscript
# Prevent empirical PSF + P finesample>1 re-entry. No fits executed.
source('R/00_load_all.R'); bd_load_all()
fixture <- matrix(1/441,21L,21L)
obs <- list(psf=fixture,psf_analytic=NULL)
spec <- list(modellist=list(pointsource=list(xcen=10,ycen=10,mag=19)),
             tofit=list(pointsource=list(xcen=TRUE,ycen=TRUE,mag=TRUE)),
             tolog=list(pointsource=list(xcen=FALSE,ycen=FALSE,mag=FALSE)))
ok<-tryCatch({bd_setup_profit_psf(spec,obs,'P',5L);FALSE},error=function(e)
  grepl('requires finesample=1',conditionMessage(e),fixed=TRUE))
stopifnot(ok)
emp <- bd_setup_profit_psf(spec,obs,'P',1L)
stopifnot(identical(emp$psf_type,'EMPIRICAL_IMAGE'),
          identical(emp$flag,'EMPIRICAL_P_FS1_RENDERER_LIMITATION'))
obs$psf_analytic<-list(fwhm=3,con=2.1,axrat=1,ang=0)
ana <- bd_setup_profit_psf(spec,obs,'P',5L)
stopifnot(identical(ana$psf_type,'ANALYTIC_MOFFAT'),
          is.null(ana$extra$psf), all(ana$extra$psfdim>=41),
          !is.null(ana$spec$modellist$psf$moffat),
          !any(unlist(ana$spec$tofit$psf)))
cat('PSF_BRANCH_CONTRACT_PASS\n')
