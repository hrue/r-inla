#!/usr/bin/env bash
## Install every system package needed to build the inla binary on
## Ubuntu/Debian. Used by CI (.github/workflows/build-inla.yml) and
## runnable as-is on any Ubuntu machine to prepare a build environment.
##
## R packages: r-base-dev provides the headers and the shared libR that
## INLA_WITH_LIBR links (rgeneric); r-mathlib is the standalone Rmath
## library used for distributions throughout. The geo libraries
## (udunits2/gdal/geos/proj) are not used by inla itself: they are runtime
## requirements of the sf/units R packages, which current fmesher imports,
## and the INLA package imports fmesher.
set -e

## No sudo when already root (e.g. inside a container).
SUDO=""
[ "$(id -u)" != 0 ] && SUDO=sudo

## Latest R from CRAN's own repository: the distro archive freezes R at
## release time (Ubuntu 24.04 carries 4.3), while users install current R
## from CRAN -- so the binary should link against current R as well.
## Ubuntu only; other apt systems keep their distro R.
if [ -r /etc/os-release ]; then . /etc/os-release; fi
if [ "${ID:-}" = "ubuntu" ] && [ ! -f /etc/apt/sources.list.d/cran.list ]; then
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | $SUDO tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc > /dev/null
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${VERSION_CODENAME}-cran40/" \
        | $SUDO tee /etc/apt/sources.list.d/cran.list > /dev/null
fi

$SUDO apt-get update
$SUDO apt-get install -y --no-install-recommends \
    gcc g++ gfortran make \
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
    libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    libnuma-dev libhwloc-dev libltdl-dev

## The pinned newer toolchain for the primary CI job. Ubuntu-archive only:
## other systems (e.g. the gcc:latest canary container) use their own
## default compiler.
if [ "${ID:-}" = "ubuntu" ]; then
    $SUDO apt-get install -y --no-install-recommends gcc-14 g++-14 gfortran-14
fi
