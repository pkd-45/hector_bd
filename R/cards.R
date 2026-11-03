# Visual diagnostic; complete card only after fitting and G3.
bd_save_card <- function(ip,fits,metrics,dir) {
  path<-file.path(dir,'card.png')
  if(file.exists(path)) return(invisible(path))
  if(!dir.exists(dir)) dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  model_name<-names(fits)[1];model<-fits[[model_name]]$model
  x<-ip$built$input$image
  stamp_region<-ip$built$input$region
  file<-tempfile(pattern='.card-',tmpdir=dir,fileext='.png')
  on.exit(unlink(file),add=TRUE)
  grDevices::png(file,width=1500,height=440,res=130)
  old<-graphics::par(mfrow=c(1,4),mar=c(3,3,3,2))
  plot_matrix<-function(z,label) {
    graphics::image(seq_len(nrow(z)),seq_len(ncol(z)),z,
      xlab='Native axis 1',ylab='Native axis 2',main=label)
    graphics::contour(seq_len(nrow(z)),seq_len(ncol(z)),
      stamp_region*1,levels=.5,drawlabels=FALSE,add=TRUE,col='red')
  }
  plot_matrix(x,'Image + original segment')
  plot_matrix(model,paste('Best',model_name))
  plot_matrix(x-ip$sky$residual_sky-model,'Residual')
  graphics::plot.new();graphics::title('Diagnostics')
  graphics::text(.05,.85,paste('Fitted',sum(stamp_region),'pixels'),adj=0)
  graphics::text(.05,.72,paste('Sky RMS',format(ip$sky$skyRMS,digits=4)),adj=0)
  graphics::text(.05,.59,paste('PSF sum',format(sum(ip$built$input$psf),digits=6)),adj=0)
  graphics::text(.05,.46,'No segment dilation or clipping',adj=0)
  graphics::par(old);grDevices::dev.off()
  if(!file.rename(file,path)) stop('Card atomic commit failed')
  invisible(path)
}
# Notebook previews use the same stage outputs; no science logic in cells.
bd_plot_stage <- function(ctx,stage) {
  if(is.null(ctx$stages[[stage]])) stop('Stage not complete: ',stage)
  if(stage %in% c('inputs','g2')) {
    ip<-ctx$stages$inputs$result
    graphics::image(ip$built$input$image, main='Original labelled support (red contour)',
                    xlab='R/FITS axis 1',ylab='R/FITS axis 2')
    graphics::contour(ip$built$input$region*1,levels=.5,drawlabels=FALSE,add=TRUE,col='red')
  } else if(stage=='f1_f9') {
    print(do.call(rbind,lapply(ctx$stages$f1_f9$result$results,as.data.frame)))
  } else if(stage=='fit') {
    fits<-ctx$stages$fit$result
    for(n in names(fits)) graphics::image(fits[[n]]$model,main=n)
  } else if(stage=='postprocess') print(ctx$stages$postprocess$result)
  else if(stage=='card') cat('Card:',file.path(ctx$outdir,'card.png'),'\n')
  invisible(TRUE)
}
