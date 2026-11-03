# Requires Rscript, package digest. No real files or fitting.
source('src/segment_contract.R')
source('src/build_inputs.R')
step0 <- list(segmentation_product='regrouped',mapping_evidence='SYNTHETIC_TEST',
 science_band='r', image_units='synthetic_counts',sky_state='synthetic_zero_sky',
 science_image_path='SYNTHETIC',segmap_path='SYNTHETIC',
 science_hdu='SYNTHETIC',segmap_hdu='SYNTHETIC',
 sigma_source='SYNTHETIC',gain_source='SYNTHETIC',wcs_verified=TRUE)
L<-matrix(0L,73,88);L[14:34,20:41]<-9L;L[56:59,65:68]<-9L;L[38:41,44:48]<-3L
img<-matrix(1,73,88);sig<-matrix(2,73,88);mask<-matrix(FALSE,73,88)
img[15,20]<-NA_real_;sig[56,66]<-0;mask[57,67]<-TRUE
L0<-L;img0<-img;sig0<-sig;mask0<-mask
z<-segment_contract(L,c(9),img,sig,mask,step0)
stopifnot(sum(z$T_g)==21*22+16,sum(z$F_g)==21*22+16-3,
          nrow(z$ledger)==3, all(z$ledger$row_y0==z$ledger$row_1based-1L))
# The remote disconnected island survives in the support.
stopifnot(z$F_g[56,65],z$F_g[59,68])
psf<-matrix(1,21,21)
b<-bd_build_inputs(L,c(9),img,sig,mask,psf,2.5,step0)
stopifnot(identical(b$input$segim,L[b$coords$row_1based[1]:b$coords$row_1based[2],
 b$coords$col_1based[1]:b$coords$col_1based[2],drop=FALSE]))
bd_assert_inputs_unchanged(b,L,img,sig,mask)
stopifnot(identical(L,L0),identical(img,img0),identical(sig,sig0),identical(mask,mask0))
fail<-try(bd_build_inputs(L,9,img,sig,mask,psf,2.5,
  within(step0,{wcs_verified<-FALSE})),silent=TRUE)
stopifnot(inherits(fail,'try-error'))
fail<-try(bd_build_inputs(L,c(9),img[-1,,drop=FALSE],sig,mask,psf,2.5,step0),silent=TRUE)
stopifnot(inherits(fail,'try-error'))
cat('D2+D3 synthetic tests PASS\n')

# Amended A1/A2: explicit delivered map subtype may be unresolved; exact
# WCS position lookup is allowed but must flag sky/boundary/segID mismatch.
step0_unclassified <- step0;step0_unclassified$segmentation_product <- 'as_delivered_unclassified'
stopifnot(isTRUE(validate_step0(step0_unclassified)))
Lw <- matrix(0L,15L,18L);Lw[4:10,5:12] <- 11L; Lw[13:15,15:18] <- 12L
wcs <- function(ra,dec) list(row=as.integer(ra),col=as.integer(dec),
  frame='ORIGINAL_SEGMAP',convention='R_1BASED_ROW_COL')
attr(wcs,'verified_original_segmap_wcs') <- TRUE
a <- bd_resolve_original_labels(Lw,7,8,wcs,catalogue_seg_id=11)
stopifnot(identical(a$A_g,11L),is.finite(a$centroid_offset_px))
ag <- bd_resolve_original_labels(Lw,7,8,wcs,catalogue_seg_id=11,
    catalogue_group_id='G1',verified_group_id='G1',group_members=c(11L,12L))
stopifnot(identical(ag$A_g,c(11L,12L)))
for (fun in list(function()bd_resolve_original_labels(Lw,2,2,wcs),
  function()bd_resolve_original_labels(Lw,4,5,wcs),
  function()bd_resolve_original_labels(Lw,7,8,wcs,catalogue_seg_id=12),
  function()bd_resolve_original_labels(Lw,7,8,wcs,catalogue_group_id='G1'))) {
   err <- try(fun(),silent=TRUE);stopifnot(inherits(err,'try-error'))
}
cat('Step0 A1/A2 synthetic rules PASS\n')
