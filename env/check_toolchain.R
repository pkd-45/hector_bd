#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
cat('R:', R.version.string, '\n')
work <- file.path(tempdir(), 'native_toolchain')
dir.create(work)
setwd(work)
writeLines('void bd_c(double *x) { *x += 1.0; }', 'probe_c.c')
writeLines('extern "C" void bd_cpp(double *x) { *x += 2.0; }', 'probe_cpp.cpp')
writeLines(c('      subroutine bd_fortran(x)', '      double precision x',
             '      x = x + 3.0d0', '      end'), 'probe_fortran.f')
status <- system2(file.path(R.home('bin'), 'R'),
                  c('CMD', 'SHLIB', 'probe_c.c', 'probe_cpp.cpp', 'probe_fortran.f',
                    '-o', 'native_probe.so'))
stopifnot(status == 0L)
dll <- dyn.load(file.path(work, 'native_probe.so'))
stopifnot(identical(.C('bd_c', x=0)$x, 1),
          identical(.C('bd_cpp', x=0)$x, 2),
          identical(.Fortran('bd_fortran', x=0)$x, 3))
dyn.unload(dll[['path']])
stopifnot(system2(Sys.getenv('AR'), c('rcs', 'native_probe.a', 'probe_c.o')) == 0L,
          system2(Sys.getenv('RANLIB'), 'native_probe.a') == 0L)
cat('NATIVE_C_CPP_FORTRAN_AR_RANLIB_VERIFIED\n')
