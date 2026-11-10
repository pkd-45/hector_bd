# End-to-end SYNTHETIC_ONLY recovery; never reads Friday/WAVES files.
# Independently reconstructed: original Claude recovery test is NOT attached.
# Random Gaussian noise; fixed seed; no deterministic stripe contamination.
source('R/00_load_all.R'); bd_load_all()
stopifnot(requireNamespace('ProFit',quietly=TRUE))
set.seed(109)
d <- c(57L,67L); yy <- 1:d[2]; xx <- 1:d[1]
X <- outer(xx,rep(1,d[2]));Y <- outer(rep(1,d[1]),yy)
truth <- list(sersic=list(xcen=c(28.7,28.7),ycen=c(33.3,33.3),
  mag=c(18.6,17.9),re=c(3.5,14),nser=c(4,1),axrat=c(1,.55),
  ang=c(0,32),box=c(0,0)))
psf <- as.matrix(utils::read.csv('tests/fixtures/f2_moffat_21px.csv',header=FALSE))
stopifnot(identical(dim(psf),c(21L,21L)),abs(sum(psf)-1)<1e-12)
# Analytic Moffat PSF centred at d/2 in native coordinates.
# No profitMakeModel PSF with a half-pixel-shifted centre.
stopifnot(which.max(psf)==which.max(.bd_f2_analytic_moffat()))
model <- ProFit::profitMakeModel(modellist=truth,dim=d,psf=psf,magzero=23.9)
if(is.list(model)) model <- model$z
stopifnot(is.matrix(model),identical(dim(model),d))
noise <- matrix(stats::rnorm(prod(d),sd=.12),d[1],d[2])
image <- model + noise
L <- matrix(0L,d[1],d[2]); L[(X-29)^2+(Y-34)^2<24^2] <- 101L
wcs <- function(ra,dec) list(row=29L,col=34L,frame='ORIGINAL_SEGMAP',
                               convention='R_1BASED_ROW_COL')
attr(wcs,'verified_original_segmap_wcs') <- TRUE
labels <- bd_resolve_original_labels(L,0,0,wcs,catalogue_seg_id=101L)
stopifnot(identical(labels$A_g,101L))
sky <- bd_sky_noise(image,L,effective_gain=687,sky_mask_radius_px=5L)
step0 <- list(segmentation_product='as_delivered_unclassified',
 mapping_evidence='SYNTHETIC_WCS',science_band='r',image_units='synthetic_microJy',
 sky_state='synthetic',science_image_path='SYNTHETIC',segmap_path='SYNTHETIC',
 science_hdu='SYNTHETIC',segmap_hdu='SYNTHETIC',sigma_source='SYNTHETIC',
 gain_source='SYNTHETIC',wcs_verified=TRUE)
mask <- matrix(FALSE,d[1],d[2]);
built <- bd_build_inputs(L,labels$A_g,image,sky$sigma,mask,psf,3,step0)
bd_assert_inputs_unchanged(built,L,image,sky$sigma,mask)
ipt <- list(built=built,observed=list(L=L,image=image,psf=psf,detector_mask=mask,
  source_hashes=list(image='SYNTHETIC')),sky=sky)
pilot <- list(models=c('S','P','D_FIXED4','BD_FREE_N'),
  magzero=23.9,finesample=1L,seeds=list(109L),optim_maxit=500L,
  synthetic_test_mode=TRUE)
start <- Sys.time()
diagnostic_dir <- Sys.getenv('BD_SYNTHETIC_REPORT_DIR',unset='')
fits <- bd_fit_four_models(ipt,pilot,diagnostic_dir=
  if(nzchar(diagnostic_dir)) file.path(diagnostic_dir,'optimizer_starts') else NULL)
elapsed <- as.numeric(difftime(Sys.time(),start,units='secs'))
stopifnot(setequal(names(fits),pilot$models),
  all(vapply(fits,function(x)x$convergence==0L && x$failed_evaluations==0L,logical(1))))
fixed <- fits$D_FIXED4$modellist_best$sersic
err_angle <- function(a,b) abs(((a-b+90) %% 180)-90)
# Reasonable recovery tolerances on a deliberately noisy finite stamp.
stopifnot(abs(fixed$xcen[1]-28.7)<1, abs(fixed$ycen[1]-33.3)<1,
 abs(fixed$mag[1]-18.6)<.65,abs(fixed$mag[2]-17.9)<.65,
 abs(fixed$re[1]-3.5)<1.5,abs(fixed$re[2]-14)<4,
 abs(fixed$axrat[2]-.55)<.18,err_angle(fixed$ang[2],32)<15)
stopifnot(fits$BD_FREE_N$objective <= fits$D_FIXED4$objective+1e-5,
          fits$D_FIXED4$objective<fits$S$objective,
          fits$D_FIXED4$objective<fits$P$objective,
          !is.null(fits$BD_FREE_N$starts$nested_from_D_FIXED4))
cat('SYNTHETIC_RECOVERY_PASS elapsed_seconds=',elapsed,
    ' D_FIXED4=',fits$D_FIXED4$objective,
    ' BD_FREE_N=',fits$BD_FREE_N$objective,'\n',sep='')
