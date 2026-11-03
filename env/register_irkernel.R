#!/usr/bin/env Rscript
# Uses the active R with the separate pinned restart library.
source(file.path('R','00_load_all.R'))
root<-bd_find_root()
.libPaths(c(file.path(root,'env','restart_Rlib'),.libPaths()))
if(!requireNamespace('IRkernel',quietly=TRUE)) stop('Install IRkernel for the active R; see the Mac setup in README.md')
IRkernel::installspec(name='waves_r_bd',displayname='R (waves_r BD pinned)',user=TRUE)
