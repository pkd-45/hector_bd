# Read verified configuration from the project; never guess unknown inputs.
bd_config <- function(profile='mac', root=bd_find_root()) {
  stopifnot(profile %in% c('mac','friday'))
  if (!requireNamespace('yaml',quietly=TRUE)) stop('R package yaml required')
  read <- function(s) yaml::read_yaml(file.path(root,'config',s))
  out <- list(root=root,profile=profile, step0=read('step0.yaml'),
              paths=read(paste0('paths_',profile,'.yaml')),
              pilot=read('run_pilot.yaml'),pins=read('package_pins.yaml'))
  out$paths$output_root <- bd_resolve_path(out$paths$output_root,root)
  out
}
bd_resolve_path <- function(path,root) {
  if (!is.character(path) || length(path)!=1L || !nzchar(path)) stop('Path missing in config')
  if (startsWith(path,'~/')) return(path.expand(path))
  if (startsWith(path,'/')) return(path)
  file.path(root,path)
}
bd_assert_output <- function(path, root, data_roots=c()) {
  # Prevent writes into historical input directories; prohibit symlink traversal.
  current <- path
  repeat {
    if (file.exists(current) && nzchar(Sys.readlink(current))) stop('Symlinked output path refused')
    parent <- dirname(current)
    if (identical(parent,current)) break
    current <- parent
  }
  here <- normalizePath(root,mustWork=TRUE)
  out <- normalizePath(dirname(path),mustWork=FALSE)
  forbidden <- unique(data_roots)
  for (f in forbidden) if (identical(out,f) || startsWith(out,paste0(f,'/'))) stop('Output in historical/input tree refused')
  if (identical(out,here) || startsWith(out,paste0(here,'/'))) {
    if (!startsWith(out,file.path(here,'outputs'))) stop('Project writes restricted to outputs/ subdir')
  }
  invisible(TRUE)
}
bd_fingerprint <- function(config) {
  packages <- names(config$pins$packages)
  tab <- lapply(packages,function(n) {
    want <- config$pins$packages[[n]]
    if (!requireNamespace(n,quietly=TRUE)) stop('Package absent: ',n)
    desc <- utils::packageDescription(n)
    remote <- desc$RemoteSha %||% ''
    if (!identical(as.character(remote),as.character(want$sha)))
      stop('Pinned SHA mismatch for ',n,': installed=',remote,' expected=',want$sha)
    if (!is.null(want$version) && as.character(utils::packageVersion(n)) != as.character(want$version))
      stop('Pinned version mismatch for ',n)
    list(name=n, version=as.character(utils::packageVersion(n)), sha=remote)
  })
  names(tab) <- packages
  tab
}
bd_source_hash <- function(root=bd_find_root()) {
  files <- system2('git',c('-C',root,'ls-files','--',
                           'R','config','env','scripts','slurm','python'),
                   stdout=TRUE,stderr=TRUE)
  status <- attr(files,'status')
  if(!is.null(status) && status!=0) stop('Unable to enumerate Git-tracked source files')
  if(!length(files)) stop('No Git-tracked source files found')
  hashes <- vapply(files,function(path) {
    digest::digest(file=file.path(root,path),algo='sha256')
  },character(1))
  digest::digest(paste(files,hashes,sep=':',collapse='\n'),algo='sha256')
}
`%||%` <- function(x,y) if(is.null(x) || length(x)==0L || is.na(x[1])) y else x
bd_file_hash <- function(path) {
  if(!file.exists(path)) stop('Missing file: ',path)
  if(!requireNamespace('digest',quietly=TRUE)) stop('digest required')
  digest::digest(file=path,algo='sha256')
}
bd_atomic_rds <- function(object, path) {
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old,object)) stop('Resume conflict; refusing to overwrite: ',path)
    return(invisible('RESUMED_IDENTICAL'))
  }
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(pattern='.bd-write-',tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(object,tmp,version=3,compress=FALSE)
  if (!file.rename(tmp,path)) stop('Atomic commit failed: ',path)
  invisible('CREATED')
}
