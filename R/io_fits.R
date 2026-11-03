# FITS and WCS adapter. Explicit HDU, original grid, no image resampling.
bd_read_fits <- function(file, hdu) {
  stopifnot(is.character(file), length(file)==1L, file.exists(file),
            length(hdu)==1L, is.numeric(hdu), hdu>=0, hdu==floor(hdu))
  if (grepl('[.]gz$',file,ignore.case=TRUE))
    stop('BLOCKED: compressed FITS decompression side effects unverified; use a verified read-only uncompressed science copy')
  if(!requireNamespace('Rfits',quietly=TRUE)) stop('Rfits missing')
  f <- Rfits::Rfits_read_image
  arguments <- names(formals(f))
  hd <- if ('hdu' %in% arguments) 'hdu' else if ('ext' %in% arguments) 'ext' else if ('extnum' %in% arguments) 'extnum' else NULL
  if (is.null(hd) && hdu!=0) stop('Explicit HDU cannot be selected with installed Rfits API')
  # FITS HDUs in the index are zero-based; Rfits ext=1 is FITS HDU 0.
  header <- Rfits::Rfits_read_header(filename=file,ext=as.integer(hdu)+1L)
  if(is.null(header$keyvalues$NAXIS) || header$keyvalues$NAXIS!=2L)
    stop('Declared HDU is not a 2D image; Rfits auto-extension fallback forbidden')
  args <- list(file); if (!is.null(hd)) args[[hd]] <- as.integer(hdu)+1L
  obj <- do.call(f,args)
  if (!is.list(obj) || !is.matrix(obj$imDat) || is.null(obj$keyvalues))
    stop('Rfits returned unexpected image/header representation; inspect installed API')
  list(z=obj$imDat, keyvalues=obj$keyvalues, file=file, hdu=hdu, sha256=bd_file_hash(file))
}
bd_check_original_grid <- function(science, segmentation) {
  if(!identical(dim(science$z),dim(segmentation$z))) stop('Original FITS grids differ')
  if(!requireNamespace('Rwcs',quietly=TRUE)) stop('Rwcs missing')
  dimxy <- dim(science$z)
  pts <- rbind(c(1,1),c(dimxy[1],1),c(1,dimxy[2]),dimxy/2,dimxy)
  w1 <- Rwcs::Rwcs_p2s(x=pts[,1],y=pts[,2],keyvalues=science$keyvalues)
  w2 <- Rwcs::Rwcs_p2s(x=pts[,1],y=pts[,2],keyvalues=segmentation$keyvalues)
  # Explicitly validate expected result columns; do not infer swapped axes.
  if (!is.matrix(w1) || !is.matrix(w2) ||
      !all(c('RA','Dec') %in% colnames(w1)) || !all(c('RA','Dec') %in% colnames(w2)))
    stop('Rwcs_p2s result conventions not verified; inspect API before running')
  dra <- ((w1[,'RA']-w2[,'RA']+180)%%360)-180
  sep <- sqrt((dra*cos(w1[,'Dec']*pi/180))^2 + (w1[,'Dec']-w2[,'Dec'])^2)*3600
  if (any(!is.finite(sep)) || max(sep)>0.01) stop('WCS disagreement over original grid >0.01 arcsec')
  invisible(list(max_wcs_delta_arcsec=max(sep),grid=dimxy))
}
bd_segmap_converter <- function(seg) {
  force(seg)
  f <- function(ra_deg,dec_deg) {
    p <- Rwcs::Rwcs_s2p(RA=ra_deg,Dec=dec_deg,keyvalues=seg$keyvalues)
    if (!is.matrix(p) || !all(c('x','y') %in% colnames(p))) stop('Unrecognised Rwcs_s2p orientation')
    x <- as.numeric(p[1,'x']);y <- as.numeric(p[1,'y'])
    # FITS pixel centers are integers; Rfits matrix axis1 corresponds to FITS x.
    if (!all(is.finite(c(x,y)))) stop('Nonfinite WCS conversion')
    list(row=as.integer(floor(x+0.5)),col=as.integer(floor(y+0.5)),
         frame='ORIGINAL_SEGMAP',convention='R_1BASED_ROW_COL')
  }
  attr(f,'verified_original_segmap_wcs') <- TRUE
  f
}
# `path_index.tsv`: one row per catalogue object, exact absolute input paths,
# generated only after Stage 0 evidence, not from plausible filenames.
bd_index_row <- function(config, id) {
  index <- bd_resolve_path(config$paths$input_index,config$root)
  if (!file.exists(index)) stop('BLOCKED: import verified Stage 0 path index: ',index)
  d <- utils::read.delim(index,stringsAsFactors=FALSE,check.names=FALSE)
  needed <- c('uberID','image_path','segmap_path','image_hdu','segmap_hdu',
              'psf_path','ra_deg','dec_deg','segID','psf_fwhm_px',
              'detector_mask_policy')
  if (!all(needed %in% names(d))) stop('Missing path index columns: ',paste(setdiff(needed,names(d)),collapse=','))
  row <- d[d$uberID==id,,drop=FALSE]
  if(nrow(row)!=1L) stop('Missing/duplicated frozen galaxy in path index: ',id)
  for (n in c('image_path','segmap_path','psf_path')) {
    if (is.na(row[[n]]) || !nzchar(row[[n]]) || !file.exists(row[[n]])) stop('Unverified/unavailable source: ',n)
    row[[paste0(n,'_sha')]] <- bd_file_hash(row[[n]])
  }
  row
}
bd_load_observed <- function(config,id) {
  row <- bd_index_row(config,id)
  sc <- bd_read_fits(row$image_path,as.integer(row$image_hdu))
  sg <- bd_read_fits(row$segmap_path,as.integer(row$segmap_hdu))
  evidence <- bd_check_original_grid(sc,sg)
  L <- sg$z
  if(any(!is.finite(L)) || any(L<0) || any(L!=floor(L))) stop('Noninteger WAVES label map')
  # No saved PSF may bypass an explicit source and kind; v14 provenance
  # must be separately signed off in the evidence registry.
  pobj <- readRDS(row$psf_path)
  if (!is.matrix(pobj)) stop('PSF RDS must be an explicitly selected matrix; inspect archive before using list-based PSFs')
  psf <- pobj
  analytic <- NULL
  kind <- attr(pobj,'kind',exact=TRUE)
  if(is.character(kind) && length(kind)==1L && grepl('moffat',kind,ignore.case=TRUE)) {
    attrs <- list(fwhm=attr(pobj,'fwhm',exact=TRUE),
                  con=attr(pobj,'beta',exact=TRUE),
                  axrat=attr(pobj,'axrat',exact=TRUE),
                  ang=attr(pobj,'ang',exact=TRUE))
    if(all(vapply(attrs,function(x)is.numeric(x)&&length(x)==1L&&is.finite(x),logical(1))))
      analytic <- attrs
    else stop('Analytic PSF kind has incomplete v14 Moffat parameters')
  }
  if (!all(is.finite(psf)) || any(psf<0) || any(dim(psf)%%2!=1)) stop('Invalid PSF')
  segid <- if (is.na(row$segID) || !nzchar(row$segID)) NULL else as.numeric(row$segID)
  groupid <- if ('groupID' %in% names(row) && !is.na(row$groupID) && nzchar(row$groupID)) row$groupID else NULL
  members <- if (!is.null(groupid)) {
    if (!'group_members_verified' %in% names(row) || row$group_members_verified!='TRUE')
      stop('Group ID membership not independently verified')
    as.integer(strsplit(row$group_members,';',fixed=TRUE)[[1]])
  } else NULL
  A <- bd_resolve_original_labels(L,as.numeric(row$ra_deg),as.numeric(row$dec_deg),
              bd_segmap_converter(sg),catalogue_seg_id=segid,
              catalogue_group_id=groupid,verified_group_id=groupid,group_members=members)
  mask <- if (identical(row$detector_mask_policy,'NONE_VERIFIED_FROM_SOURCE'))
      matrix(FALSE,nrow(L),ncol(L))
    else stop('Detector mask policy unresolved; no default blank mask')
  fw <- as.numeric(row$psf_fwhm_px)
  if(!is.finite(fw) || fw<=0) stop('Per-object measured PSF FWHM missing')
  # Exact source SHA plus WCS check supplies G0 evidence when authenticated.
  list(L=L,image=sc$z,psf=psf,psf_analytic=analytic,labels=A,source_files=row,
       grid_evidence=evidence, detector_mask=mask,psf_fwhm=fw,
       source_hashes=list(image=sc$sha256,segmap=sg$sha256,psf=bd_file_hash(row$psf_path)))
}
