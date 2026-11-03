#!/usr/bin/env Rscript
# Synthetic, four-quadrant ProFit angle-convention regression; no observed data.
source('R/00_load_all.R');bd_load_all()
stopifnot(requireNamespace('ProFit',quietly=TRUE))
d<-c(151L,151L); centre<-75.5
reg<-matrix(TRUE,d[1],d[2]); err<-numeric()
angle_diff<-function(a,b) abs(((a-b+90) %% 180)-90)
for(pa in c(-80,-44,-14,18,46,79,112,149)) {
  ml<-list(sersic=list(xcen=centre,ycen=centre,mag=17,re=17,nser=1,
                        axrat=.43,ang=pa,box=0))
  z<-ProFit::profitMakeModel(modellist=ml,dim=d,magzero=25,psf=NULL)
  if(is.list(z))z<-z$z
  stopifnot(is.matrix(z),all(is.finite(z)))
  measured<-bd_segment_moments(z,reg)
  e<-angle_diff(measured$pa,pa)
  cat('MOMENT_ANGLE',pa,'measured',measured$pa,'error_deg',e,'\n')
  err<-c(err,e)
}
stopifnot(all(err<2))
cat('MOMENT_PA_4_QUADRANTS_PASS max_error_deg=',max(err),'\n',sep='')
