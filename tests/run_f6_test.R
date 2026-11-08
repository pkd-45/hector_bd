#!/usr/bin/env Rscript
# Isolated synthetic F6 diagnostic. Run through SLURM after install_pinned.R.
if (!nzchar(Sys.getenv('SLURM_JOB_ID')))
  stop('F6 rendering must run in a SLURM allocation')
source(file.path('R', '00_load_all.R'))
root <- bd_find_root()
lib <- normalizePath(file.path(root, 'env', 'restart_Rlib'), mustWork=TRUE)
.libPaths(c(lib, .libPaths()))
if (!requireNamespace('yaml', quietly=TRUE)) stop('yaml is required to read package pins')
pins <- yaml::read_yaml(file.path(root, 'config', 'package_pins.yaml'))$packages
for (package in c('ProFit', 'digest')) {
  wanted <- pins[[package]]
  desc <- suppressWarnings(utils::packageDescription(package, lib.loc=lib))
  if (!is.list(desc) || !identical(desc$RemoteSha, wanted$sha))
    stop('F6 requires the exact pinned source SHA for ', package)
  if (!is.null(wanted$version) &&
      !identical(as.character(utils::packageVersion(package, lib.loc=lib)),
                 as.character(wanted$version)))
    stop('F6 requires the pinned version for ', package)
  if (!requireNamespace(package, quietly=TRUE)) stop('Cannot load ', package)
  path <- normalizePath(getNamespaceInfo(asNamespace(package), 'path'), mustWork=TRUE)
  if (!identical(path, normalizePath(file.path(lib, package), mustWork=TRUE)))
    stop('F6 loaded ', package, ' from outside restart_Rlib: ', path)
  cat('F6_PACKAGE', package, 'VERSION', desc$Version, 'REMOTE_SHA', desc$RemoteSha,
      'LOADED_PATH', path, '\n')
}
cat('F6_RUNTIME', R.version.string, 'R_HOME', R.home(),
    'SLURM_JOB_ID', Sys.getenv('SLURM_JOB_ID'), '\n')
cat('F6_LIB_PATHS', paste(.libPaths(), collapse=';'), '\n')
dlls <- getNamespaceInfo(asNamespace('ProFit'), 'DLLs')
if (!length(dlls)) stop('No loaded ProFit shared library found')
for (dll in dlls) {
  path <- normalizePath(dll[['path']], mustWork=TRUE)
  if (!startsWith(path, paste0(lib, '/ProFit/')))
    stop('Loaded ProFit shared library is outside restart_Rlib: ', path)
  cat('F6_PROFIT_BINARY_SHA256', path,
      digest::digest(file=path, algo='sha256'), '\n')
}
sys.source(file.path(root, 'R', 'renderer_tests.R'), envir=.GlobalEnv)
d <- c(87L, 109L); xc <- 27.75; yc <- 70.25
psf <- ProFit::profitMakeGaussianPSF(fwhm=3, dim=c(21L, 21L))
if (is.list(psf) && !is.null(psf$z)) psf <- psf$z
psf <- psf/sum(psf)
# Keep an independent baseline for comparison; pass the original PSF unchanged
# to every renderer call, so this cannot conceal a native input-mutation bug.
baseline <- unserialize(serialize(psf, NULL))
state <- function(x) list(
  object_sha256=digest::digest(x, algo='sha256'),
  values_sha256=digest::digest(writeBin(as.vector(x), raw(), size=8L,
                                       endian='little'),
                               algo='sha256', serialize=FALSE),
  attributes=attributes(x), sum=sum(x), range=range(x))
initial <- state(psf)
cat('F6_INITIAL_PSF_SHA256', initial$object_sha256, '\n')
cat('F6_INITIAL_VALUES_SHA256', initial$values_sha256, '\n')
cat('F6_INITIAL_ATTRIBUTES\n'); dput(initial$attributes)
first_change <- NA_integer_
first_value_change <- NA_integer_
first_attribute_change <- NA_integer_
after_first <- NULL
maximum_difference <- 0
for (i in seq_len(100L)) {
  invisible(.bd_model(d, xc, yc, psf=psf))
  current <- state(psf)
  difference <- max(abs(psf-baseline))
  maximum_difference <- max(maximum_difference, difference)
  if (is.na(first_change) &&
      !identical(initial$object_sha256, current$object_sha256)) first_change <- i
  if (is.na(first_value_change) &&
      !identical(initial$values_sha256, current$values_sha256)) first_value_change <- i
  if (is.na(first_attribute_change) &&
      !identical(initial$attributes, current$attributes)) first_attribute_change <- i
  if (i == 1L) {
    after_first <- current
    cat('F6_AFTER_FIRST_PSF_SHA256', current$object_sha256, '\n')
    cat('F6_AFTER_FIRST_VALUES_SHA256', current$values_sha256, '\n')
    cat('F6_AFTER_FIRST_MAX_ABS_DIFF', sprintf('%.17g', difference), '\n')
    cat('F6_AFTER_FIRST_ATTRIBUTES\n'); dput(current$attributes)
  }
}
final <- state(psf)
cat('F6_FINAL_PSF_SHA256', final$object_sha256, '\n')
cat('F6_FINAL_VALUES_SHA256', final$values_sha256, '\n')
cat('F6_FINAL_MAX_ABS_DIFF', sprintf('%.17g', max(abs(psf-baseline))), '\n')
cat('F6_MAX_ABS_DIFF_ALL_CALLS', sprintf('%.17g', maximum_difference), '\n')
cat('F6_FIRST_CHANGED_RENDER', first_change, '\n')
cat('F6_FIRST_VALUES_CHANGED_RENDER', first_value_change, '\n')
cat('F6_FIRST_ATTRIBUTES_CHANGED_RENDER', first_attribute_change, '\n')
cat('F6_FINAL_ATTRIBUTES\n'); dput(final$attributes)
passed <- identical(initial$object_sha256, final$object_sha256)
cat('F6', if (passed) 'PASS' else 'FAIL', '\n')
if (!passed) stop('F6 input PSF SHA-256 changed after profitMakeModel calls')
