# Measured segment-radius photometry with explicit finite-stamp status.
bd_postprocess <- function(ip,fits) {
  rr<-ip$built$input$region
  if (!is.matrix(rr)) stop('G2 region missing')
  pos<-which(rr,arr.ind=TRUE); cx<-mean(pos[,1]);cy<-mean(pos[,2])
  n<-nrow(pos); radius=sqrt(n/pi)
  out<-lapply(names(fits),function(name) {
    model<-fits[[name]]$intrinsic_model
    if(!identical(dim(model),dim(rr))) stop('Model/native stamp mismatch')
    flux_segment<-sum(model[rr]);
    x<-outer(seq_len(nrow(rr)),rep(1,ncol(rr))); y<-outer(rep(1,nrow(rr)),seq_len(ncol(rr)))
    aperture <- (x-cx)^2 + (y-cy)^2 <= radius^2
    touches_edge <- any(aperture[c(1,nrow(aperture)),]) || any(aperture[,c(1,ncol(aperture))])
    flux_radius <- if(touches_edge) NA_real_ else sum(model[aperture])
    data.frame(model=name,segment_pixels=n,segment_radius_equiv_px=radius,
       intrinsic_flux_inside_original_label=flux_segment,
       intrinsic_flux_inside_equiv_radius=flux_radius,
       aperture_clipped=touches_edge,
       size_definition='EQUIVALENT_CIRCULAR_RADIUS; NOT SERSIC SEMI_MAJOR_SEGRAD',
       stringsAsFactors=FALSE)
  })
  do.call(rbind,out)
}
bd_collate <- function(output_root,run_id) {
  here<-file.path(output_root,run_id)
  if(!dir.exists(here)) stop('No run directory')
  entries<-list.dirs(here,recursive=FALSE,full.names=TRUE)
  rows<-lapply(entries,function(d) {
    path<-file.path(d,'postprocess.rds')
    if (!file.exists(path)) return(NULL)
    z<-readRDS(path)$result
    z$uberID<-basename(d);z
  })
  if (!length(Filter(Negate(is.null),rows))) stop('No complete results to collate')
  out<-do.call(rbind,Filter(Negate(is.null),rows))
  output<-file.path(here,'collated.tsv')
  if(file.exists(output)) stop('Collation exists; no overwrite')
  utils::write.table(out,output,sep='\t',row.names=FALSE,quote=TRUE)
  output
}
