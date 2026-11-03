source(file.path('R','00_load_all.R'));bd_load_all()
L<-matrix(0L,25,25);L[5:15,5:15]<-10L
f<-function(ra,dec)list(row=as.integer(ra),col=as.integer(dec),frame='ORIGINAL_SEGMAP',convention='R_1BASED_ROW_COL')
attr(f,'verified_original_segmap_wcs')<-TRUE
normal<-bd_resolve_original_labels(L,10,10,f,catalogue_seg_id=10)
stopifnot(identical(normal$A_g,10L),identical(normal$boundary_flag,FALSE))
edge<-bd_resolve_original_labels(L,5,5,f,catalogue_seg_id=10)
stopifnot(identical(edge$A_g,10L),identical(edge$boundary_flag,TRUE))
for(g in list(function()bd_resolve_original_labels(L,2,2,f),
              function()bd_resolve_original_labels(L,1,1,f),
              function()bd_resolve_original_labels(L,10,10,f,catalogue_seg_id=11))) {
  stopifnot(inherits(try(g(),silent=TRUE),'try-error'))
}
cat('AMENDMENT A1 label boundary FLAGS rather than STOPS: PASS\n')
