# D6: F1-F9 synthetic regression harness. NO REAL FITTING.
# Loaded by R/00_load_all.R
# Runs only in the pinned R/ProFit environment; a failure or an unrun test
# cannot be promoted into G4 PASS. The result defaults to BLOCKED.
.bd_test <- function(name, fn) {
  tryCatch({r <- fn();
    if (!isTRUE(r)) stop('test returned non-TRUE');
    list(gate=name,status='PASS',detail='')},
    error=function(e) list(gate=name,status='FAIL',detail=conditionMessage(e)))
}
.bd_model <- function(dim,xc,yc,re=9,nser=1,q=.55,ang=37,psf=NULL,finesample=1L) {
  m <- list(sersic=list(xcen=xc,ycen=yc,mag=17,re=re,nser=nser,
                        axrat=q,ang=ang,box=0))
  z <- ProFit::profitMakeModel(modellist=m,dim=as.integer(dim),psf=psf,
                               magzero=25,finesample=as.integer(finesample))
  if (is.list(z) && !is.null(z$z)) z$z else z
}
.bd_ll <- function(image,model,sigma,region) {
  # Exact norm objective used in current ProFit: dnorm(standardised residual).
  # This is cross-checked against profitLikeModel below, NOT assumed equivalent.
  sum(stats::dnorm(((image-model)/sigma)[region], log=TRUE))
}
.bd_pro_like <- function(image, model, sigma, region, psf, ml) {
  dat <- ProFit::profitSetupData(image=image,region=region,sigma=sigma,
    segim=matrix(0L,nrow(image),ncol(image)),mask=matrix(FALSE,nrow(image),ncol(image)),
    modellist=ml,tofit=list(sersic=list(xcen=FALSE,ycen=FALSE,mag=FALSE,
      re=FALSE,nser=FALSE,axrat=FALSE,ang=FALSE,box=FALSE)),
    psf=psf,magzero=25,like.func='norm',algo.func='LA',finesample=1L,
    nbenchmark=0L)
  stopifnot(identical(dat$region,region))
  ans <- ProFit::profitLikeModel(parm=numeric(0),Data=dat,model=list(z=model))
  -ans$Dev/2
}
run_renderer_tests <- function() {
  stopifnot(requireNamespace('ProFit',quietly=TRUE),requireNamespace('digest',quietly=TRUE))
  d <- c(87L,109L);xc <- 27.75;yc <- 70.25
  # This PSF is already normalised; use same object on all tests.
  psf <- ProFit::profitMakeGaussianPSF(fwhm=3,dim=c(21L,21L))
  if (is.list(psf) && !is.null(psf$z)) psf <- psf$z
  psf <- psf/sum(psf)
  region <- matrix(FALSE,d[1],d[2]);region[9:54,35:96] <- TRUE
  model <- .bd_model(d,xc,yc,psf=psf)
  image <- model + .1; sigma <- matrix(2,d[1],d[2])
  ml <- list(sersic=list(xcen=xc,ycen=yc,mag=17,re=9,nser=1,
      axrat=.55,ang=37,box=0))
  list(
    F1=.bd_test('F1',function(){
      pad<-32L;big<-.bd_model(d+2*pad,xc+pad,yc+pad,psf=psf)
      z<-big[(pad+1):(pad+d[1]),(pad+1):(pad+d[2])]
      # Compare on fitted support in noise units, not full image edges.
      max(abs((z-model)[region])/sigma[region])<.01
    }),
    F2=.bd_test('F2',function(){ bd_f2_run() }),
    F3=.bd_test('F3',function(){
      z0<-.bd_model(d,xc,yc,re=7,q=.33,ang=0,psf=NULL)
      z9<-.bd_model(d,xc,yc,re=7,q=.33,ang=90,psf=NULL)
      x <- seq_len(d[1])-.5; y <- seq_len(d[2])-.5
      v1<-function(z) sum(z*outer((x-xc)^2,rep(1,d[2])))/sum(z)
      v2<-function(z) sum(z*outer(rep(1,d[1]),(y-yc)^2))/sum(z)
      v2(z0)>2*v1(z0) && v1(z9)>2*v2(z9)
    }),
    F4=.bd_test('F4',function(){
      z<-.bd_model(d,xc,yc,re=5,ang=0,q=1,psf=NULL)
      mr<-sum(z*outer(seq_len(d[1])-.5,rep(1,d[2])))/sum(z)
      mc<-sum(z*outer(rep(1,d[1]),seq_len(d[2])-.5))/sum(z)
      abs(mr-xc)<.1 && abs(mc-yc)<.1
    }),
    F5=.bd_test('F5',function(){
      # Pixel-integrated ProFit n=1 Sérsic: elliptical aperture with semi-major
      # re should contain approximately half the total large-stamp flux.
      re<-16; q<-.7; bx<-200.5;by<-200.5; dim1<-c(401L,401L)
      z<-.bd_model(dim1,bx,by,re=re,nser=1,q=q,ang=0,psf=NULL)
      xx<-seq_len(dim1[1])-.5;yy<-seq_len(dim1[2])-.5
      r2<-outer(((xx-bx)/q)^2,rep(1,length(yy)))+
          outer(rep(1,length(xx)),(yy-by)^2)
      frac<-sum(z[r2<=re^2])/sum(z)
      abs(frac-.5)<.025 && abs(re/1.67835*1.67835-re)<1e-12
    }),
    F6=.bd_test('F6',function(){
      before<-digest::digest(psf,algo='sha256');
      for (i in seq_len(100)) invisible(.bd_model(d,xc,yc,psf=psf))
      identical(before,digest::digest(psf,algo='sha256'))
    }),
    F7=.bd_test('F7',function(){
      z1<-sum((((image-model)/sigma)[region])^2)
      z2<-sum((((image-model)/(2*sigma))[region])^2)
      abs(z2/z1-.25)<1e-12
    }),
    F8=.bd_test('F8',function(){
      expected<-.bd_ll(image,model,sigma,region)
      observed<-.bd_pro_like(image,model,sigma,region,psf,ml)
      is.finite(observed) && abs(expected-observed)<1e-8
    }),
    F9=.bd_test('F9',function(){
      before<-.bd_pro_like(image,model,sigma,region,psf,ml)
      outside<-image;outside[1,1]<-outside[1,1]+1000
      inside<-image;inside[20,40]<-inside[20,40]+1000
      a<-.bd_pro_like(outside,model,sigma,region,psf,ml)
      b<-.bd_pro_like(inside,model,sigma,region,psf,ml)
      identical(a,before) && b!=before
    })
  )
}
if (sys.nframe()==0L) {
  results <- run_renderer_tests()
  for (nm in names(results)) cat(nm,results[[nm]]$status,results[[nm]]$detail,'\n')
  if (!all(vapply(results,function(x) identical(x$status,'PASS'),logical(1))))
    stop('G4 BLOCKED: not every F1-F9 renderer test passed')
}
