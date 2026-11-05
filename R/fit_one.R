# Shared notebook/CLI pipeline. All stages call exactly bd_stage().
# Real fitting REQUIRES explicit approval; synthetic tests never fit objects.
bd_synthetic_galaxy <- function(id) {
  n<-65L;L<-matrix(0L,n,n)
  x<-outer(seq_len(n),rep(1L,n));y<-outer(rep(1L,n),seq_len(n))
  L[((x-33)^2+(y-33)^2)<250]<-101L
  L[8:9,58:59]<-101L  # detached island retained in G2
  image<-exp(-((x-33)^2+(y-33)^2)/45)*250 + 0.02*sin(x*y)
  psf<-outer(-10:10,-10:10,function(a,b)exp(-(a^2+b^2)/(2*1.5^2)))
  psf<-psf/sum(psf)
  list(L=L,image=image,psf=psf,labels=list(A_g=101L,boundary_flag=FALSE,
     boundary_check='3x3_UNIFORM_LABEL'),
     detector_mask=matrix(FALSE,n,n),source_hashes=list(image='SYNTHETIC',segmap='SYNTHETIC',psf='SYNTHETIC'))
}
bd_context <- function(id, config, run_id='test_synthetic', mode='synthetic', output_override=NULL) {
  if (!grepl('^W[0-9]+$',id)) stop('Invalid immutable uberID')
  if (!grepl('^[A-Za-z0-9_-]+$',run_id)) stop('Invalid immutable run_id')
  if (!mode %in% c('synthetic','observed')) stop('Unknown mode')
  outroot <- if (is.null(output_override)) config$paths$output_root else output_override
  path <- file.path(outroot,run_id,id)
  bd_assert_output(path,config$root,config$paths$data_roots %||% character(0))
  e<-new.env(parent=emptyenv());e$id<-id;e$config<-config;e$run_id<-run_id
  e$mode<-mode;e$outdir<-path;e$stages<-list();e
}
bd_step0_gate <- function(config) {
  required<-c('science_image','original_segmentation','wcs_alignment','units_zeropoint',
              'upstream_sky','gain_calibration','noise_protocol','catalogue_label_mapping')
  ev<-config$step0$evidence
  if (!all(required %in% names(ev)) ||
      !all(vapply(ev[required],function(x) identical(x$status,'VERIFIED_FROM_BYTES') &&
            nzchar(x$sha256 %||% ''),logical(1))))
    stop('G0 BLOCKED: Step 0 evidence not VERIFIED_FROM_BYTES with SHA256 for all required items')
  if(!identical(config$step0$effective_gain$status,'VERIFIED_FROM_BYTES'))
    stop('G0 BLOCKED: effective gain still only REPORTED')
  invisible(TRUE)
}
bd_psf_pilot_gate <- function(ctx) {
  # A passing synthetic F2 fixture is NOT a per-object PSF validation.
  p <- bd_resolve_path(ctx$config$pilot$per_object_psf_gate_index,ctx$config$root)
  if(!file.exists(p)) stop('F2/G5 BLOCKED: per-object PSF gate manifest absent')
  tab <- utils::read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
  fields <- c('uberID','psf_sha256','f2_status','e1_e5_status','evidence_log','evidence_sha256')
  if(!all(fields %in% names(tab))) stop('PSF gate manifest columns missing')
  row <- tab[tab$uberID==ctx$id,,drop=FALSE]
  if(nrow(row)!=1L || !identical(row$f2_status,'PASS_VERIFIED') ||
      !identical(row$e1_e5_status,'PASS_VERIFIED'))
    stop('Per-object F2 or PSF E1-E5 unverified for ',ctx$id)
  for(family in c('S','P','D_FIXED4','BD_FREE_N')) {
    col <- paste0('f2_',family)
    wanted <- as.integer(ctx$config$pilot$finesample_by_model[[family]] %||%
                         ctx$config$pilot$finesample)
    if(!col %in% names(row) || !is.finite(wanted) ||
       is.na(row[[col]]) || as.integer(row[[col]])!=wanted)
      stop('F2 per-object fine-sampling evidence missing or disagrees for ',
           ctx$id,' ',family)
  }

  actual <- ctx$stages$inputs$result$observed$source_hashes$psf
  if(!identical(row$psf_sha256,actual) || !file.exists(row$evidence_log) ||
     !identical(bd_file_hash(row$evidence_log),row$evidence_sha256))
    stop('Per-object PSF log/hash mismatch')
  invisible(TRUE)
}
bd_stage_input_hash <- function(res,ctx) {
  observed <- if(!is.null(res$result$observed)) res$result$observed else
    ctx$stages$inputs$result$observed
  hashes <- observed$source_hashes
  files <- observed$source_files
  if(identical(ctx$mode,'observed')) {
    required_hashes <- c('image','segmap','psf')
    if (!is.list(hashes) ||
        !identical(sort(names(hashes)),sort(required_hashes)) ||
        !all(vapply(hashes,function(x)
          is.character(x) && length(x)==1L &&
          !is.na(x) && grepl('^[0-9a-fA-F]{64}$',x),
          logical(1))))
      stop('Observed source hashes missing, malformed or incomplete')
    # Fail closed, and bind the cache to the CURRENT verified index row: bd_index_row re-reads
    # labels, coordinates, HDUs and paths and re-hashes the image, segmap and PSF files.
    if(!(is.data.frame(files) && nrow(files)==1L))
      stop('Observed stage lacks its verified index row; provenance cannot be checked; do not resume')
    current <- bd_index_row(ctx$config,ctx$id)
    if(!identical(.bd_hash(as.list(current)),.bd_hash(as.list(files))))
      stop('Verified index row changed (labels, coordinates, HDUs, paths or file bytes); choose new run_id')
  }
  if(is.data.frame(files) && nrow(files)==1L) {
    for (kind in names(hashes)) {
      path_col <- paste0(kind,'_path')
      if(!path_col %in% names(files)) stop('Input path missing for ',kind)
      current <- bd_file_hash(files[[path_col]][[1]])
      if(!identical(current,hashes[[kind]]))
        stop('Input bytes changed for ',kind)
      hashes[[kind]] <- current
    }
  }
  .bd_hash(list(hashes=hashes,index_row=if(is.data.frame(files)) as.list(files) else NULL))
}
bd_stage <- function(ctx, stage) {
  stopifnot(is.environment(ctx),stage %in% c('inputs','g2','f1_f9','fit','postprocess','card'))
  deps<-list(inputs=NULL,g2='inputs',f1_f9='g2',fit='f1_f9',postprocess='fit',card='postprocess')
  dep<-deps[[stage]]
  if (!is.null(dep) && is.null(ctx$stages[[dep]])) stop('Run stage ',dep,' before ',stage)
  artifact<-file.path(ctx$outdir,paste0(stage,'.rds'))
  if (file.exists(artifact)) {
    hashfile<-paste0(artifact,'.sha256')
    if(!file.exists(hashfile) || !identical(trimws(readLines(hashfile,warn=FALSE)[1]),bd_file_hash(artifact)))
      stop('Saved stage artifact SHA256 missing/mismatch; do not resume')
    res<-readRDS(artifact)
    if (!identical(res$stage,stage) || !identical(res$id,ctx$id) ||
        !identical(res$mode,ctx$mode) || !identical(res$config_sha,.bd_hash(ctx$config[c('step0','paths','pilot','pins')])))
       stop('Existing stage artifact conflicts with current config; choose new run_id')
    if(!identical(res$source_sha,bd_source_hash(ctx$config$root)))
      stop('Existing stage artifact was produced by different source code; choose new run_id')
    if(!identical(res$input_sha,bd_stage_input_hash(res,ctx)))
      stop('Existing stage artifact input bytes changed; choose new run_id')
    ctx$stages[[stage]]<-res
    return(res)
  }
  config<-ctx$config
  if(stage=='inputs') {
    if(ctx$mode=='observed') {
      bd_step0_gate(config)
      obs<-bd_load_observed(config,ctx$id)
      if (isTRUE(obs$labels$boundary_flag)) warning('A1 boundary flagged; exact centre label retained')
      gain<-config$step0$effective_gain$value
      if (!is.numeric(gain) || !is.finite(gain) || gain<=0) stop('Verified effective gain required')
    } else {obs<-bd_synthetic_galaxy(ctx$id);obs$psf_fwhm<-3.5;gain<-687}
    sky<-bd_sky_noise(obs$image,obs$L,gain,
      sky_mask_radius_px=as.integer(config$step0$sky_mask_radius_px %||% 10L))
    step0 <- list(segmentation_product=config$step0$product,mapping_evidence='ORIGINAL_WCS',
      science_band='r',image_units='microJy',sky_state='measured_background',
      science_image_path=if(ctx$mode=='synthetic') 'SYNTHETIC' else obs$source_files$image_path,
      segmap_path=if(ctx$mode=='synthetic') 'SYNTHETIC' else obs$source_files$segmap_path,
      science_hdu=if(ctx$mode=='synthetic') 'SYNTHETIC' else as.character(obs$source_files$image_hdu),
      segmap_hdu=if(ctx$mode=='synthetic') 'SYNTHETIC' else as.character(obs$source_files$segmap_hdu),
      sigma_source='SKY_FLOOR_AND_MODEL',
      gain_source='VERIFIED_EFFECTIVE_GAIN',wcs_verified=TRUE)
    built<-bd_build_inputs(obs$L,as.integer(obs$labels$A_g),obs$image,
      sky$sigma,obs$detector_mask,obs$psf,fwhm_px=obs$psf_fwhm,step0=step0)
    out<-list(observed=obs,sky=sky,built=built)
  } else if(stage=='g2') {
    ip<-ctx$stages$inputs$result
    bd_assert_inputs_unchanged(ip$built,ip$observed$L,ip$observed$image,
                               ip$sky$sigma,ip$observed$detector_mask)
    out<-list(gate='G2_PASS',full_coords=ip$built$contract$hashes$support_coords,
              excluded=nrow(ip$built$contract$ledger),labels=ip$observed$labels)
  } else if(stage=='f1_f9') {
    # Gate must execute on the pinned R runtime; failing even one test blocks fit.
    fingerprint<-bd_fingerprint(config)
    tests<-run_renderer_tests()
    passed<-all(vapply(tests,function(x)identical(x$status,'PASS'),logical(1)))
    if (!passed) stop('G4 BLOCKED: F1-F9 failed; no fit allowed')
    out<-list(gate='G4_PASS',results=tests,packages=fingerprint)
  } else if(stage=='fit') {
    if(ctx$mode!='observed') stop('Synthetic mode stops before fitting')
    if (!identical(config$pilot$pilot_approved,TRUE) ||
        !identical(config$pilot$psf_g5_approved,TRUE) ||
        !identical(config$pilot$models_g6_approved,TRUE))
      stop('Fitting BLOCKED: approval, PSF G5 and model G6 must be explicitly TRUE')
    if (!ctx$id %in% config$pilot$galaxy_ids) stop('Galaxy not in frozen 20-object pilot')
    if(length(config$pilot$galaxy_ids)!=20L || anyDuplicated(config$pilot$galaxy_ids))
      stop('Pilot must be a frozen list of exactly 20 unique real galaxy IDs')
    bd_psf_pilot_gate(ctx)
    out<-bd_fit_four_models(ctx$stages$inputs$result,config$pilot,
       verified_gate_token=.bd_hash(list(ctx$id,ctx$stages$g2$result,
           ctx$stages$f1_f9$result,config$pilot$pilot_approved)),
       diagnostic_dir=file.path(ctx$outdir,'optim_diagnostics'))
  } else if(stage=='postprocess') {
    out<-bd_postprocess(ctx$stages$inputs$result,ctx$stages$fit$result)
  } else if(stage=='card') {
    bd_save_card(ctx$stages$inputs$result,ctx$stages$fit$result,
                 ctx$stages$postprocess$result,ctx$outdir)
    out<-list(path='card.png')
  }
  res<-list(stage=stage,id=ctx$id,mode=ctx$mode,
            config_sha=.bd_hash(config[c('step0','paths','pilot','pins')]),
            source_sha=bd_source_hash(config$root),
            input_sha=bd_stage_input_hash(list(result=out),ctx),
            result=out)
  bd_atomic_rds(res,artifact)
  sha<-bd_file_hash(artifact)
  hashfile<-paste0(artifact,'.sha256')
  if(file.exists(hashfile)) {
    if(!identical(trimws(readLines(hashfile,warn=FALSE)[1]),sha)) stop('Stage checksum conflicts')
  } else writeLines(sha,hashfile,useBytes=TRUE)
  ctx$stages[[stage]]<-res
  res
}
bd_run_one <- function(id, config,run_id='test_synthetic',mode='synthetic',
                       output_override=NULL,stages=c('inputs','g2','f1_f9')) {
  ctx<-bd_context(id,config,run_id,mode,output_override)
  for (stage in stages) bd_stage(ctx,stage)
  ctx
}
# The exact same native ProFit evaluation path serves the notebook and CLI.
bd_render_profit <- function(parm, dat) {
  check <- dat
  check$algo.func <- 'check'
  rendered <- ProFit::profitLikeModel(parm=parm,Data=check)
  if(!is.list(rendered) || !is.list(rendered$model) ||
     !is.matrix(rendered$model$z) ||
     !identical(dim(rendered$model$z),dim(dat$image)) ||
     any(!is.finite(rendered$model$z)))
    stop('ProFit returned invalid model image; refusing objective evaluation')
  rendered$model$z
}
bd_nested_free_n_ml <- function(d_fixed_ml, bd_free_spec) {
  x <- bd_free_spec$modellist
  stopifnot(!is.null(d_fixed_ml$sersic),length(d_fixed_ml$sersic$mag)==2L)
  for(field in names(x$sersic)) x$sersic[[field]] <- d_fixed_ml$sersic[[field]]
  # Explicit exact embedding of D_FIXED4 in BD_FREE_N.
  x$sersic$nser <- c(4,1)
  x$sersic$axrat[1] <- 1
  x$sersic$ang[1] <- 0
  x
}
# Analytic PSF is the PSF-v14 Moffat *parameter model*, not a resampled
# image masquerading as an analytic kernel. All values must be from verified
# provenance (RDS attributes for saved v14 PSFs or synthetic fixture).
bd_psf_model <- function(observed) {
  a <- observed$psf_analytic
  if (is.null(a)) return(NULL)
  required <- c('fwhm','con','axrat','ang')
  if (!is.list(a) || !all(required %in% names(a)) ||
      any(!is.finite(as.numeric(unlist(a[required])))) ||
      a$fwhm<=0 || a$con<=1 || a$axrat<=0 || a$axrat>1)
    stop('Invalid analytic Moffat PSF provenance')
  list(moffat=list(mag=0,fwhm=a$fwhm,con=a$con,
                   axrat=a$axrat,ang=a$ang,box=0))
}
bd_setup_profit_psf <- function(spec, observed, family, fs) {
  analytic <- bd_psf_model(observed)
  if (!is.null(analytic)) {
    spec$modellist$psf <- analytic
    spec$tofit$psf <- .bd_flags(analytic)
    spec$tolog$psf <- .bd_flags(analytic)
    # Need a broad PSF domain for the ProFit internal flux-normalisation guard.
    # psfdim is native pixels, not a pre-fine-sampled bitmap.
    d <- max(81L,as.integer(max(dim(observed$psf))));if(d%%2L==0L)d<-d+1L
    return(list(spec=spec,extra=list(psf=NULL,psfdim=c(d,d)),
      psf_type='ANALYTIC_MOFFAT',finesample=fs,flag=''))
  }
  # P at fs>1 is broken in pinned ProFit for empirical bitmap PSFs.
  if (family=='P' && fs!=1L)
    stop('P empirical PSF requires finesample=1; cannot use bitmap with finesample>1')
  list(spec=spec,extra=list(psf=observed$psf),
       psf_type='EMPIRICAL_IMAGE',finesample=fs,
       flag=if(family=='P') 'EMPIRICAL_P_FS1_RENDERER_LIMITATION' else '')
}
bd_fit_four_models <- function(ip, pilot, verified_gate_token=NULL,diagnostic_dir=NULL) {
  if(isTRUE(pilot$synthetic_test_mode)) {
    if(!identical(ip$observed$source_hashes$image,'SYNTHETIC'))
      stop('Synthetic fit test must not read observed image data')
  } else if(!is.character(verified_gate_token) ||
            length(verified_gate_token)!=1L || nchar(verified_gate_token)!=64L)
    stop('Observed fitter cannot be called without stage-gate provenance')
  # Callable only after manually checked G0/G2/G4/G5/G6 (in bd_stage).
  built <- ip$built
  bd_assert_inputs_unchanged(built,ip$observed$L,ip$observed$image,
                              ip$sky$sigma,ip$observed$detector_mask)
  z <- built$input$image-ip$sky$residual_sky
  region <- built$input$region
  if (!identical(region,built$contract$F_g[
       built$coords$row_1based[1]:built$coords$row_1based[2],
       built$coords$col_1based[1]:built$coords$col_1based[2],drop=FALSE]))
    stop('G2 drift before fit')
  required <- c('S','P','D_FIXED4','BD_FREE_N')
  if(!setequal(pilot$models,required) || anyDuplicated(pilot$models))
    stop('Exactly the four frozen model families are required')
  # The nesting direction is scientifically required; no order-only bug.
  families <- required
  output <- list()
  for (family in families) {
    spec <- bd_model_spec(family,built,pilot$magzero)
    family_fs <- pilot$finesample_by_model[[family]] %||% pilot$finesample
    if(!is.numeric(family_fs) || length(family_fs)!=1L ||
       !family_fs %in% c(1L,3L,5L,7L,9L)) stop('Unverified family fine-sampling')
    psf_obs <- ip$observed; psf_obs$psf <- built$input$psf   # use the normalised, hash-checked PSF
    psf_mode <- bd_setup_profit_psf(spec,psf_obs,family,as.integer(family_fs))
    spec <- psf_mode$spec
    dat <- do.call(ProFit::profitSetupData,c(list(image=z,region=region,sigma=built$input$sigma,
       segim=built$input$segim,mask=built$input$mask,
       modellist=spec$modellist,tofit=spec$tofit,tolog=spec$tolog,
       intervals=NULL,constraints=spec$constraints,
       psffinesampled=FALSE,finesample=as.integer(family_fs),
       magzero=pilot$magzero,like.func='norm',algo.func='optim',nbenchmark=0L,
       verbose=FALSE),psf_mode$extra))
    if (!identical(dat$psftype,if(psf_mode$psf_type=='ANALYTIC_MOFFAT') 'analytical' else 'empirical'))
      stop('ProFit PSF type differs from explicit source classification')
    if(!identical(dat$region,region)) stop('G2: ProFit changed region')
    if(length(dat$init)!=spec$k) stop('Parameter vector mismatch for ',family)
    base <- bd_extract_free(spec$modellist,spec$tofit)
    if(!isTRUE(all.equal(as.numeric(dat$init),base,tolerance=1e-10)))
      stop('ProFit free parameter ordering/transform does not match declared model')
    state <- new.env(parent=emptyenv());state$evaluations<-0L;state$failed<-0L
    objective <- function(p) {
      state$evaluations <- state$evaluations + 1L
      mu <- tryCatch(bd_render_profit(p,dat),error=function(e) {
        state$failed <- state$failed + 1L
        stop(sprintf('%s render failed at evaluation %d: %s',family,
                     state$evaluations,conditionMessage(e)),call.=FALSE)
      })
      sigma <- bd_model_sigma(mu,ip$sky$skyRMS,ip$sky$gain)
      if(any(!is.finite(sigma[region])) || any(sigma[region]<=0)) {
        state$failed <- state$failed + 1L
        stop('Invalid model-dependent sigma in objective: ',family)
      }
      result <- sum(.5*((z[region]-mu[region])/sigma[region])^2 +
                    log(sigma[region]) + .5*log(2*pi))
      if(!is.finite(result)) {
        state$failed <- state$failed + 1L
        stop('Non-finite likelihood for ',family)
      }
      result
    }
    starts <- bd_structured_modellists(spec$modellist,family)
    if(isTRUE(pilot$synthetic_test_mode))
      starts <- starts[c('base','re_half','bt_bulge_bright')]
    if(family=='BD_FREE_N') {
      if(is.null(output$D_FIXED4) || is.null(output$D_FIXED4$modellist_best))
        stop('Nested warm start missing converged D_FIXED4 model')
      starts$nested_from_D_FIXED4 <- bd_nested_free_n_ml(
        output$D_FIXED4$modellist_best,spec)
    }
    records <- list()
    start_hashes <- character()
    for(start_name in names(starts)) {
      p <- bd_extract_free(starts[[start_name]],spec$tofit)
      if(length(p)!=spec$k || any(!is.finite(p)))
        stop('Malformed structured start for ',family,': ',start_name)
      p <- pmin(pmax(p,spec$bounds$lower),spec$bounds$upper)
      h <- .bd_hash(p)
      if(h %in% start_hashes) next
      start_hashes <- c(start_hashes,h)
      # Preflight render. Never turn a failed render into a constant objective plateau.
      initial_value <- objective(p)
      fit <- stats::optim(p,objective,method='L-BFGS-B',
        lower=spec$bounds$lower,upper=spec$bounds$upper,
        control=list(maxit=as.integer(pilot$optim_maxit)))
      fit$start <- p;fit$start_objective <- initial_value
      fit$start_name <- start_name
      fit$failed_evaluations <- state$failed
      fit$evaluations_total <- state$evaluations
      fit$stationary <- isTRUE(all.equal(as.numeric(fit$par),p,tolerance=1e-10))
      if(!is.finite(fit$value) || state$failed>0)
        stop('Unresolved objective failure during ',family,' ',start_name)
      records[[start_name]] <- fit
    }
    # Persist every start, including line-search (code 52) and stationary failures.
    # An uncertified result is diagnostic evidence, NEVER a science fit.
    uncertified <- lapply(Filter(function(x) x$convergence!=0L || isTRUE(x$stationary),records),
      function(v) list(start_name=v$start_name,parameters=v$par,
        objective=v$value,convergence=v$convergence,
        message=v$message %||% '',bound_hits=which(
          abs(v$par-spec$bounds$lower)<1e-10 |
          abs(v$par-spec$bounds$upper)<1e-10),
        stationary=v$stationary))
    if (!is.null(diagnostic_dir)) {
      if(!dir.exists(diagnostic_dir)) dir.create(diagnostic_dir,recursive=TRUE)
      # Never silently drop failed starts or promote them to a certified fit.
      dpath <- file.path(diagnostic_dir,paste0(family,'_all_starts.rds'))
      bd_atomic_rds(list(family=family,all=records,uncertified=uncertified,
                   science_result='NOT_CERTIFIED_BY_THIS_LOG'),dpath)
    }
    converged <- Filter(function(x) x$convergence==0L && !isTRUE(x$stationary),records)

    if(length(converged)==0L) stop('No converged non-stationary fit for ',family)
    ord <- order(vapply(converged,function(x)x$value,numeric(1)))
    best <- converged[[ord[1]]]
    if(isTRUE(best$stationary)) stop('Best fit did not move from its start: ',family)
    # BD_FREE_N must not underperform the nested D_FIXED4 solution.
    if(family=='BD_FREE_N') {
      nested <- records[['nested_from_D_FIXED4']]
      if(is.null(nested)) stop('Nested warm-start candidate missing')
      if(nested$start_objective > output$D_FIXED4$objective + 1e-5)
        stop('D_FIXED4 cannot be embedded at identical likelihood in BD_FREE_N')
      if(best$value > output$D_FIXED4$objective + 1e-5)
        stop('G6: BD_FREE_N worse than nested D_FIXED4, local minimum unresolved')
    }
    mod <- bd_render_profit(best$par,dat)
    best_ml <- ProFit::profitRemakeModellist(parm=best$par,Data=dat)$modellist
    raw_ml <- best_ml; raw_ml$psf <- NULL
    rawmodel <- ProFit::profitMakeModel(modellist=raw_ml,dim=as.integer(dim(z)),
                                       magzero=pilot$magzero,psf=NULL)
    if(is.list(rawmodel) && !is.null(rawmodel$z)) rawmodel <- rawmodel$z
    if(!is.matrix(rawmodel) || !identical(dim(rawmodel),dim(z)))
      stop('Intrinsic light profile dimensions invalid')
    output[[family]] <- list(parameters=best$par,objective=best$value,
      model=mod,intrinsic_model=rawmodel,modellist_best=best_ml,
      k=spec$k,model_spec=spec$modelhash,
      chosen_start=best$start_name,finesample=family_fs,
      psf_type=psf_mode$psf_type,psf_renderer_flag=psf_mode$flag,
      convergence=best$convergence,convergence_message=best$message %||% '',
      evaluations_total=state$evaluations,failed_evaluations=state$failed,
      uncertified_starts=uncertified,
      starts=lapply(records,function(v)list(name=v$start_name,
        initial_objective=v$start_objective,objective=v$value,
        convergence=v$convergence,message=v$message %||% '',
        stationary=v$stationary,evaluations=v$counts)),
      status='OPTIM_ONLY_NO_POSTERIOR_ERRORS',
      psf_unchanged=identical(.bd_hash(built$input$psf),built$hashes$psf))
    if(!output[[family]]$psf_unchanged) stop('G3 PSF mutated during fit')
  }
  bd_assert_inputs_unchanged(built,ip$observed$L,ip$observed$image,
                             ip$sky$sigma,ip$observed$detector_mask)
  output
}
bd_run_pilot <- function(config, run_id=config$pilot$run_id, mode='observed') {
  ids<-config$pilot$galaxy_ids
  if (length(ids)!=20L || anyDuplicated(ids)) stop('20 declared pilot IDs must be frozen before any fit')
  for (id in ids) bd_run_one(id,config,run_id,mode,
              stages=c('inputs','g2','f1_f9','fit','postprocess','card'))
  bd_collate(config$paths$output_root,run_id)
}
