#!/usr/bin/env bash
## Install every system package needed to build the inla binary on Fedora.
## Used by the newest-toolchain CI canary (Fedora tracks the newest GCC and
## current R), and runnable as-is on any Fedora machine.
set -e

SUDO=""
[ "$(id -u)" != 0 ] && SUDO=sudo

$SUDO dnf -y install \
    gcc gcc-c++ gcc-gfortran make git-core findutils diffutils rsync \
    openblas-devel \
    gsl-devel \
    metis-devel \
    muParser-devel \
    openssl-devel \
    suitesparse-devel \
    zlib-devel \
    eigen3-devel \
    simde-devel \
    R-core-devel libRmath-devel \
    udunits2-devel gdal-devel geos-devel proj-devel \
    numactl-devel hwloc-devel libtool-ltdl-devel
