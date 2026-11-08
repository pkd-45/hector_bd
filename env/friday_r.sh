#!/usr/bin/env bash
# Source from a Friday SLURM job. Reuse the existing native toolchain without
# conda activate (its R hook runs javareconf in the historical environment).
BD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
BD_CONDA_PREFIX=/home/uqpdas/.conda/envs/waves_r
export PATH="$BD_CONDA_PREFIX/bin:$PATH"
export R_LIBS_USER="$BD_ROOT/env/restart_Rlib"
export R_LIBS_SITE="$BD_CONDA_PREFIX/lib/R/library"
unset R_LIBS
# Ignore inherited temporary system-compiler overrides and user startup files.
export R_MAKEVARS_USER=/dev/null
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export CC=x86_64-conda-linux-gnu-cc
export CXX=x86_64-conda-linux-gnu-c++
export FC=x86_64-conda-linux-gnu-gfortran
export F77="$FC"
export AR=x86_64-conda-linux-gnu-ar
export RANLIB=x86_64-conda-linux-gnu-ranlib
export NM=x86_64-conda-linux-gnu-nm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-1}"
for BD_TOOL in "$CC" "$CXX" "$FC" "$AR" "$RANLIB" "$NM"; do
  [[ "$(command -v "$BD_TOOL")" == "$BD_CONDA_PREFIX/bin/$BD_TOOL" ]] || {
    echo "Missing native Conda tool: $BD_TOOL" >&2
    return 1
  }
done
unset BD_TOOL
