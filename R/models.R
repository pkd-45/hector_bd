# Frozen four-family definitions; all coordinates in ProFit native axes.
.bd_flags <- function(x) lapply(x, function(comp) lapply(comp, function(v) rep(FALSE, length(v))))
.bd_bounds <- function(ml, lo, hi, tofit) {
  # Flattening preserves ProFit's modellist component + parameter order.
  keep <- unlist(tofit,use.names=FALSE)
  list(lower=as.numeric(unlist(lo,use.names=FALSE)[keep]),
       upper=as.numeric(unlist(hi,use.names=FALSE)[keep]))
}
# Robust positive-flux, within-original-F_g second moments. Returned angle uses
# the ProFit convention: 0 degrees extends along the second R array axis.
bd_segment_moments <- function(image, region) {
  coords <- which(region, arr.ind=TRUE)
  stopifnot(nrow(coords)>=5L)
  values <- image[region]
  floor <- stats::median(values[is.finite(values)])
  w <- pmax(values-floor, 0)
  if (sum(w)<=0 || !all(is.finite(w))) w <- rep(1,length(values))
  xy <- cbind(coords[,1]-.5, coords[,2]-.5)
  centre <- colSums(xy*w)/sum(w)
  delta <- sweep(xy,2,centre)
  cov <- crossprod(delta,delta*w)/sum(w)
  eig <- eigen(cov,symmetric=TRUE)
  q <- max(.12,min(1,sqrt(max(eig$values[2],1e-7)/max(eig$values[1],1e-7))))
  pa <- -atan2(eig$vectors[1,1],eig$vectors[2,1])*180/pi
  pa <- ((pa+90) %% 180)-90
  list(xcen=centre[1],ycen=centre[2],q=q,pa=pa,
       covariance=cov,flux=sum(w),n_pixels=nrow(coords),
       convention='ProFit_0_MAJOR_AXIS_SECOND_R_DIM')
}
# Free parameters in model field/component order; fail if ProFit disagrees.
bd_extract_free <- function(ml,tofit) {
  stopifnot(identical(names(ml),names(tofit)))
  out <- numeric()
  for (comp in names(ml)) for (field in names(ml[[comp]])) {
    flag <- tofit[[comp]][[field]]
    stopifnot(length(flag)==length(ml[[comp]][[field]]))
    out <- c(out,as.numeric(ml[[comp]][[field]][flag]))
  }
  out
}
# Structured variants, never small uncalibrated Gaussian parameter noise.
bd_structured_modellists <- function(ml, family) {
  tweak <- function(name,fn){ x<-ml; x<-fn(x); x }
  shift_pa <- function(x, delta){
    if (!is.null(x$sersic)) x$sersic$ang <- ((x$sersic$ang+delta+180) %% 360)-180
    x
  }
  scale_re <- function(x, factor){
    if (!is.null(x$sersic)) x$sersic$re <- x$sersic$re*factor
    x
  }
  flux_split <- function(x, bulge_bright){
    if (family %in% c('D_FIXED4','BD_FREE_N')) {
      # MAG shifts preserve order, and deliberately explore B/T degeneracy.
      x$sersic$mag <- x$sersic$mag + if(bulge_bright)c(-.75,.75) else c(.75,-.75)
    } else if(family=='P') {
      x$pointsource$mag <- x$pointsource$mag + if(bulge_bright)-.75 else .75
      x$sersic$mag <- x$sersic$mag + if(bulge_bright).75 else -.75
    }
    x
  }
  list(base=ml,
       pa_minus45=shift_pa(ml,-45),pa_plus45=shift_pa(ml,45),
       re_half=scale_re(ml,.5),re_double=scale_re(ml,2),
       bt_bulge_bright=flux_split(ml,TRUE),
       bt_disc_bright=flux_split(ml,FALSE))
}
bd_model_spec <- function(family,built,magzero) {
  stopifnot(family %in% c('S','P','D_FIXED4','BD_FREE_N'))
  region <- built$input$region
  pos <- which(region,arr.ind=TRUE); z<-built$input$image
  geom <- bd_segment_moments(z,region)
  cx <- geom$xcen; cy <- geom$ycen
  re <- max(1.5,sqrt(nrow(pos)/pi)/2)
  flux <- max(1e-4,sum(pmax(z[region],0)))
  mag <- magzero-2.5*log10(flux)
  one <- function(m, radius, n, q=.75, pa=0) list(xcen=cx,ycen=cy,mag=m,re=radius,
                                           nser=n,axrat=q,ang=pa,box=0)
  if (family=='S') ml<-list(sersic=one(mag,re,2,geom$q,geom$pa))
  if (family=='P') ml<-list(pointsource=list(xcen=cx,ycen=cy,mag=mag+1),
                       sersic=one(mag+.55,re,1,geom$q,geom$pa))
  if (family %in% c('D_FIXED4','BD_FREE_N')) {
    b<-one(mag+1,re*.4,if(family=='D_FIXED4') 4 else 2,
           if(family=='D_FIXED4') 1 else geom$q,
           if(family=='D_FIXED4') 0 else geom$pa)
    disc<-one(mag+.6,re*1.6,1,geom$q,geom$pa)
    ml<-list(sersic=list(xcen=c(cx,cx),ycen=c(cy,cy),mag=c(mag+1,mag+.6),
      re=c(re*.4,re*1.6),nser=c(if(family=='D_FIXED4')4 else 2,1),
      axrat=c(if(family=='D_FIXED4')1 else geom$q,geom$q),
      ang=c(if(family=='D_FIXED4')0 else geom$pa,geom$pa),box=c(0,0)))
  }
  tf<-.bd_flags(ml)
  if (family=='S') for(k in c('xcen','ycen','mag','re','nser','axrat','ang')) tf$sersic[[k]]<-TRUE
  if (family=='P') {
    for(k in names(tf$pointsource)) tf$pointsource[[k]]<-TRUE
    for(k in c('mag','re','axrat','ang')) tf$sersic[[k]]<-TRUE
  }
  if (family=='D_FIXED4') {
    for(k in c('xcen','ycen')) tf$sersic[[k]]<-c(TRUE,FALSE)
    for(k in c('mag','re')) tf$sersic[[k]]<-c(TRUE,TRUE)
    for(k in c('axrat','ang')) tf$sersic[[k]]<-c(FALSE,TRUE)
  }
  if (family=='BD_FREE_N') {
    for(k in c('xcen','ycen')) tf$sersic[[k]]<-c(TRUE,FALSE)
    for(k in c('mag','re','axrat','ang')) tf$sersic[[k]]<-c(TRUE,TRUE)
    tf$sersic$nser<-c(TRUE,FALSE)
  }
  tie <- switch(family,
    P=function(m){m$sersic$xcen<-m$pointsource$xcen;m$sersic$ycen<-m$pointsource$ycen;m},
    D_FIXED4=function(m){m$sersic$xcen[2]<-m$sersic$xcen[1];m$sersic$ycen[2]<-m$sersic$ycen[1];m},
    BD_FREE_N=function(m){m$sersic$xcen[2]<-m$sersic$xcen[1];m$sersic$ycen[2]<-m$sersic$ycen[1];m},
    NULL)
  lo <- lapply(ml,function(x)lapply(x,function(v)rep(-Inf,length(v))))
  hi <- lapply(ml,function(x)lapply(x,function(v)rep(Inf,length(v))))
  for(comp in names(ml)) {
    for (field in names(ml[[comp]])) {
      if(field=='xcen'){lo[[comp]][[field]][]<-0;hi[[comp]][[field]][]<-nrow(z)}
      if(field=='ycen'){lo[[comp]][[field]][]<-0;hi[[comp]][[field]][]<-ncol(z)}
      if(field=='mag'){lo[[comp]][[field]][]<-mag-8;hi[[comp]][[field]][]<-mag+8}
      if(field=='re'){lo[[comp]][[field]][]<-.3;hi[[comp]][[field]][]<-max(dim(z))}
      if(field=='nser'){lo[[comp]][[field]][]<-.3;hi[[comp]][[field]][]<-8}
      if(field=='axrat'){lo[[comp]][[field]][]<-.1;hi[[comp]][[field]][]<-1}
      if(field=='ang'){lo[[comp]][[field]][]<--180;hi[[comp]][[field]][]<-180}
    }
  }
  bounds <- .bd_bounds(ml,lo,hi,tf)
  k<-length(bounds$lower)
  expected<-c(S=7L,P=7L,D_FIXED4=8L,BD_FREE_N=11L)[[family]]
  if(k!=expected) stop('Incorrect declared free parameter count in ',family)
  list(family=family,modellist=ml,tofit=tf,tolog=.bd_flags(ml),
       constraints=tie,bounds=bounds,k=k,priors='NONE',
       geom=geom,
       semantics='D_FIXED4 round n4 bulge; P tied centres; BD_FREE_N tied centres',
       modelhash=.bd_hash(list(ml=ml,tofit=tf,bounds=bounds,constraints=deparse(tie))))
}
