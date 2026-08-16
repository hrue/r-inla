#!/usr/bin/env bash
## Install every system package needed to build the inla binary on
## Ubuntu/Debian. Used by CI (.github/workflows/build-inla.yml) and
## runnable as-is on any Ubuntu machine to prepare a build environment.
##
## R packages: r-base-dev provides the headers and the shared libR that
## INLA_WITH_LIBR links (rgeneric); r-mathlib is the standalone Rmath
## library used for distributions throughout.
set -e

## No sudo when already root (e.g. inside a container).
SUDO=""
[ "$(id -u)" != 0 ] && SUDO=sudo

$SUDO apt-get update
$SUDO apt-get install -y --no-install-recommends \
    gcc-14 g++-14 gfortran-14 gcc g++ gfortran make \
    libopenblas-dev \
    libgsl-dev \
    libmetis-dev \
    libmuparser-dev \
    libssl-dev \
    libsuitesparse-dev \
    zlib1g-dev \
    libeigen3-dev \
    libsimde-dev \
    r-base-dev r-mathlib \
    libnuma-dev libhwloc-dev libltdl-dev
