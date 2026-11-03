# The single source-loading entry point for notebooks, CLI, and tests.
bd_find_root <- function(start = getwd()) {
  cur <- normalizePath(start, mustWork=TRUE)
  repeat {
    if (file.exists(file.path(cur, '.hector_bd_root'))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop('No .hector_bd_root found in ancestor directories')
    cur <- parent
  }
}
bd_load_all <- function(root=bd_find_root()) {
  files <- c('config_io.R','io_fits.R','segment_contract.R','build_inputs.R',
    'sky_noise.R','psf_v14.R','psf_compare.R','models.R',
    'f2_finesample_test.R','renderer_tests.R','postprocess.R','cards.R','fit_one.R')
  for (f in files) sys.source(file.path(root, 'R', f), envir=.GlobalEnv)
  invisible(files)
}
