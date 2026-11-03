# F2: empirical extended + analytic-Moffat point-source tests on SYNTHETIC data.
# DO NOT use this harness to fit observed data; never call optim/MCMC.
# It calls the actual pinned ProFit profitSetupData and profitLikeModel renderer.
# All comparisons use identical native-grid original-region pixels and native sigma.
# F2 does NOT certify per-object measured PSFs; repeat with each pilot's adopted PSF.

.bd_f2_require <- function(cond, message) {
  if (!isTRUE(cond)) stop(message, call. = FALSE)
}
.bd_f2_ml <- function(id, re_override=NULL, ps_mag_override=NULL) {
  if (id == 'compact_n4') {
    list(sersic=list(xcen=30.1,ycen=32.45,mag=19.0,
      re=if(is.null(re_override)) 1.5 else re_override,
      nser=4,axrat=.77,ang=23,box=0))
  } else if (id == 'point_plus_exponential') {
    list(pointsource=list(xcen=30.1,ycen=32.45,
      mag=if(is.null(ps_mag_override)) 19.7 else ps_mag_override),
      sersic=list(xcen=30.1,ycen=32.45,mag=18.9,re=7.0,
        nser=1,axrat=.72,ang=32,box=0))
  } else stop('Unknown synthetic F2 model')
}
.bd_f2_fixed_tofit <- function(ml) {
  tf <- lapply(ml, function(m) lapply(m,function(v) rep(FALSE,length(v))))
  comp <- if('pointsource' %in% names(tf)) 'pointsource' else 'sersic'
  tf[[comp]]$mag[1] <- TRUE  # upstream inheritIDs zero-free-parameter bug
  tf
}
.bd_f2_render <- function(ml, finesample, psf, d, region, sig) {
  .bd_f2_require(requireNamespace('ProFit',quietly=TRUE), 'BLOCKED: ProFit unavailable')
  .bd_f2_require(finesample %in% c(1L,3L,5L,7L,9L), 'F2 restricted to 1, 3, 5, 7, 9')
  .bd_f2_require(is.matrix(psf) && all(dim(psf)==c(21L,21L)) &&
    all(is.finite(psf)) && min(psf)>=0 && abs(sum(psf)-1)<1e-12,
    'F2 pinned native PSF must be a normalized 21 x 21 matrix')
  zeros <- matrix(0, nrow=d[1], ncol=d[2])
  is_point <- 'pointsource' %in% names(ml)
  if (is_point) {
    # The empirical bitmap point-source path has a pinned ProFit rendering bug
    # for finesample>1. Use an analytic PSF model, never a bitmap interpolant.
    beta <- 2.1; alpha <- 2.23873756
    fwhm <- 2*alpha*sqrt(2^(1/beta)-1)
    ml$psf <- list(moffat=list(mag=0,fwhm=fwhm,con=beta,
                                axrat=1,ang=0,box=0))
    tf <- .bd_f2_fixed_tofit(ml)
    tf$psf <- .bd_flags(ml$psf)
    extra <- list(psf=NULL,psfdim=c(81L,81L))
  } else {
    tf <- .bd_f2_fixed_tofit(ml)
    extra <- list(psf=psf)
  }
  data <- do.call(ProFit::profitSetupData,c(list(image=zeros, region=region,
    sigma=sig,segim=matrix(1L,nrow=d[1],ncol=d[2]),
    mask=matrix(FALSE,nrow=d[1],ncol=d[2]),
    modellist=ml,tofit=tf,psffinesampled=FALSE,
    finesample=as.integer(finesample),magzero=25,like.func='norm',
    algo.func='LA',nbenchmark=0L,benchintmethods='brute',
    benchconvmethods='brute',verbose=FALSE),extra))
  .bd_f2_require(identical(data$psftype,if(is_point)'analytical' else 'empirical'),
     'F2 did not use the declared analytic/empirical renderer branch')
  .bd_f2_require(identical(as.integer(data$finesample),as.integer(finesample)),
    'ProFit changed finesample')
  .bd_f2_require(identical(as.logical(data$region),as.logical(region)),
    'ProFit changed F2 test region')
  .bd_f2_require(is.matrix(data$psf) && all(dim(data$psf)%%2L==1L),
    'ProFit PSF was not an odd-sized matrix')
  .bd_f2_require(abs(sum(data$psf)-1)<1e-9, 'ProFit PSF lost normalization')
  if (finesample>1L) .bd_f2_require(all(dim(data$psf)>dim(psf)),
    'F2 PSF was not fine-sampled by profitSetupData')
  # The check path returns the model image
  # after the same ProFit convolution machinery used by the likelihood.
  # If the installed pinned version lacks this path, FAIL closed.
  check_data <- data; check_data$algo.func <- 'check'
  out <- ProFit::profitLikeModel(parm=data$init,Data=check_data)
  .bd_f2_require(is.list(out) && is.list(out$model) && is.matrix(out$model$z),
    'Pinned ProFit did not return model through algo.func=check')
  z <- out$model$z
  .bd_f2_require(identical(dim(z),as.integer(d)) && all(is.finite(z)),
    'F2 model dimensions/values invalid')
  z
}
.bd_f2_data <- function(a,b,sig,region) {
  z <- (a[region]-b[region])/sig[region]
  list(rms=sqrt(mean(z*z)),max=max(abs(z)),dchisq=sum(z*z))
}
.bd_f2_proxy <- function(case, fs, psf, d, region, sig, reference) {
  if (case=='compact_n4') {
    values <- 1.5 + seq(-2,2)*.05
    ml_for <- function(value) .bd_f2_ml(case,re_override=value)
  } else {
    values <- 19.7 + seq(-2,2)*.01
    ml_for <- function(value) .bd_f2_ml(case,ps_mag_override=value)
  }
  scores <- vapply(values,function(value){
    z <- .bd_f2_render(ml_for(value),fs,psf,d,region,sig)
    sum(((z[region]-reference[region])/sig[region])^2)
  },numeric(1))
  list(best=values[which.min(scores)],minchi=scores[which.min(scores)],
    edge=which.min(scores)%in%c(1L,length(values)))
}
# Analytic Moffat surface-brightness PSF, evaluated on the native grid.
# This validates the fixture independently of the ProFit renderer.
.bd_f2_analytic_moffat <- function(dim=21L,alpha=2.23873756,beta=2.1) {
  center <- dim/2
  axis <- seq_len(dim)-.5-center
  r2 <- outer(axis^2,rep(1,dim))+outer(rep(1,dim),axis^2)
  psf <- (1+r2/alpha^2)^(-beta)
  psf/sum(psf)
}
# Family-specific accuracy against a high fine-sampling reference using an
# independently analytic native Moffat PSF (for extended-family fixture). This is a synthetic accuracy test,
# NOT an analytic integral of a Sérsic over each native pixel.
bd_f2_accuracy <- function(case,psf,d,region,sig) {
  reference <- .bd_f2_render(.bd_f2_ml(case),9L,psf,d,region,sig)
  metric <- lapply(c(1L,3L,5L,7L),function(fs) {
    x <- .bd_f2_render(.bd_f2_ml(case),fs,psf,d,region,sig)
    .bd_f2_data(x,reference,sig,region)
  })
  names(metric) <- c('1','3','5','7')
  for(fs in names(metric)) cat('F2_ANALYTIC_MOFFAT_REFERENCE',case,'fs',fs,
      sprintf('max_sigma=%.9g delta_chi2=%.9g',
              metric[[fs]]$max,metric[[fs]]$dchisq),'\n')
  # F2 3->5 must still pass on its own. Only FS choices checked against
  # high-res reference are eligible, with the unchanged significance budget.
  candidates <- names(metric)[vapply(metric,function(m)
    is.finite(m$max)&&m$max<.1&&m$dchisq<1,logical(1))]
  chosen <- if(length(candidates)) as.integer(candidates[[1]]) else NA_integer_
  cat('F2_ACCURACY_CHOSEN',case,if(is.na(chosen))'BLOCKED' else chosen,'\n')
  list(chosen=chosen,metrics=metric,reference=reference,
       accuracy_pass=!is.na(chosen))
}
bd_f2_run <- function() {
  .bd_f2_require(requireNamespace('digest',quietly=TRUE), 'BLOCKED: digest unavailable')
  path <- file.path('tests','fixtures','f2_moffat_21px.csv')
  .bd_f2_require(file.exists(path), 'BLOCKED: pinned F2 PSF fixture not found')
  filehash <- digest::digest(file=path,algo='sha256')
  expectedhash <- 'd569969a59f40745c323250f888e3e122a457f9eb63e78029fd684982b2eb1d0'
  .bd_f2_require(identical(filehash,expectedhash),'BLOCKED: F2 pinned PSF hash differs')
  psf <- as.matrix(utils::read.csv(path,header=FALSE,check.names=FALSE))
  original_psf_hash <- digest::digest(psf,algo='sha256')
  .bd_f2_require(identical(dim(psf),c(21L,21L)), 'F2 PSF fixture shape altered')
  .bd_f2_require(abs(sum(psf)-1)<1e-12,'F2 PSF fixture not normalized')
  analytic <- .bd_f2_analytic_moffat()
  .bd_f2_require(max(abs(psf-analytic))<1e-7,
     'Pinned Moffat fixture is incompatible with independent analytic PSF')
  d<-c(63L,67L)
  xi<-seq_len(d[1])-.5; yi<-seq_len(d[2])-.5
  dist2 <- outer((xi-30.1)^2,rep(1,d[2]))+
           outer(rep(1,d[1]),(yi-32.45)^2)
  region <- dist2<=18^2
  .bd_f2_require(sum(region)>500, 'F2 region invalid')
  # Compact case fixed noise; point-source case predeclared peak SNR=60 below.
  sig_compact <- matrix(2.5,d[1],d[2])
  allpass <- logical(0)
  pkgdir <- system.file(package='ProFit')
  .bd_f2_require(nzchar(pkgdir), 'F2 cannot locate installed ProFit files')
  rlibs <- paste(.libPaths(),collapse=';')
  cat('F2_RUNTIME',R.version.string,'ProFit',
      as.character(utils::packageVersion('ProFit')),
      'PROFIT_DIR',normalizePath(pkgdir,mustWork=TRUE),
      'LIB_PATHS',rlibs,'PSF_SHA256',filehash,'\n')
  libs <- list.files(file.path(pkgdir,'libs'),
     pattern='\\.(so|dylib|dll)$',full.names=TRUE,recursive=TRUE)
  .bd_f2_require(length(libs)>0L, 'F2 cannot fingerprint pinned ProFit shared library')
  for (lib in libs) cat('F2_PROFIT_BINARY_SHA256',
      normalizePath(lib,mustWork=TRUE),digest::digest(file=lib,algo='sha256'),'\n')
  for (case in c('compact_n4','point_plus_exponential')) {
    base <- .bd_f2_ml(case)
    if (case=='point_plus_exponential') {
      # Choose noise from the analytic-path model peak (not a tuned gate value).
      pilot_sigma <- matrix(1,d[1],d[2])
      truth_snr <- .bd_f2_render(base,9L,psf,d,region,pilot_sigma)
      sig <- matrix(max(truth_snr[region])/60,d[1],d[2])
      .bd_f2_require(is.finite(sig[1,1]) && sig[1,1]>0,'F2 SNR calibration invalid')
      cat('F2_POINT_NUCLEAR_SNR 60 SIGMA',sig[1,1],
          'PSF_MODE analytic_moffat_model','\n')
      # An otherwise uninformative half-magnitude proxy invalidates F2.
      perturbed <- .bd_f2_ml(case,ps_mag_override=19.7+.5)
      zhalf <- .bd_f2_render(perturbed,9L,psf,d,region,sig)
      identifiability <- sum(((truth_snr[region]-zhalf[region])/sig[region])^2)
      .bd_f2_require(identifiability>1,'F2 nuclear magnitude is unidentifiable at declared SNR')
    } else sig <- sig_compact
    z <- lapply(c(1L,3L,5L),function(fs) .bd_f2_render(base,fs,psf,d,region,sig))
    names(z)<-c('1','3','5')
    m13<-.bd_f2_data(z[['1']],z[['3']],sig,region)
    m35<-.bd_f2_data(z[['3']],z[['5']],sig,region)
    m15<-.bd_f2_data(z[['1']],z[['5']],sig,region)
    # Fixed small one-parameter profile scan against finesample=5 synthetic
    # truth. This is a convergence proxy, not parameter recovery from real data.
    params <- lapply(c(1L,3L,5L),function(fs)
      .bd_f2_proxy(case,fs,psf,d,region,sig,z[['5']]))
    names(params)<-c('1','3','5')
    step <- if(case=='compact_n4') .05 else .01
    param_delta <- abs(params[['3']]$best-params[['5']]$best)
    ok <- is.finite(m35$max) && m35$max<.1 &&
      m35$dchisq<1 && param_delta<=step+1e-10 &&
      !params[['3']]$edge && !params[['5']]$edge
    for (comparison in c('1to3','3to5','1to5')) {
      metrics <- switch(comparison,'1to3'=m13,'3to5'=m35,'1to5'=m15)
      cat('F2_METRIC',case,comparison,
        sprintf('rms_sigma=%.9g max_sigma=%.9g delta_chi2=%.9g',
                metrics$rms,metrics$max,metrics$dchisq),'\n')
    }
    cat('F2_PARAMETER_PROXY',case,
      sprintf('best_fs1=%.8g best_fs3=%.8g best_fs5=%.8g fs3_minus_fs5=%.8g pass=%s',
              params[['1']]$best,params[['3']]$best,params[['5']]$best,
              param_delta,if(ok)'TRUE' else 'FALSE'),'\n')
    accuracy <- bd_f2_accuracy(case,psf,d,region,sig)
    # Original 3->5 gate is never relaxed: if it fails, finesample=3 is
    # explicitly unsafe. The next (5->9) convergence gate uses the identical
    # 0.1 sigma / delta chi2 <1 thresholds against the high-res reference.
    selected <- accuracy$chosen
    if(!is.na(selected)) {
      chosen_profile <- .bd_f2_proxy(case,selected,psf,d,region,sig,
                                    accuracy$reference)
      truth_profile <- .bd_f2_proxy(case,9L,psf,d,region,sig,
                                   accuracy$reference)
      accuracy_param_pass <- abs(chosen_profile$best-truth_profile$best) <= step+1e-10 &&
          !chosen_profile$edge && !truth_profile$edge
    } else accuracy_param_pass <- FALSE
    # Compact n=4 is preregistered as a finesample 3 case and must pass 3->5;
    # the unresolved P family may instead require finesample 5, but only when
    # the independently anchored 5->9 accuracy gate passes.
    permitted <- if(case=='compact_n4') ok && !is.na(selected) && selected<=3L &&
      accuracy_param_pass else ok && accuracy$accuracy_pass && accuracy_param_pass
    cat('F2_FAMILY_DECISION',case,'FS3_CONVERGED',ok,
      'SELECTED',if(is.na(selected))'BLOCKED' else selected,
      'ACCURACY_PARAMETER_PASS',accuracy_param_pass,
      'FAMILY_PASS',permitted,'\n')
    allpass <- c(allpass,permitted)
  }
  .bd_f2_require(identical(original_psf_hash,digest::digest(psf,algo='sha256')),
      'F2 pinned native PSF mutated by renderer')
  all(allpass)
}
