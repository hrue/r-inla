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
if [ "${ID:-}" = "ubuntu" ]; then
    ## The key is installed UNCONDITIONALLY: the runner image pre-adds the CRAN
    ## repository itself, so guarding the key behind "cran.list absent" skips
    ## it exactly when the image ships the repo without its key -- and apt
    ## update then dies with NO_PUBKEY before anything else runs. Re-writing
    ## the key when it exists is a no-op.
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | $SUDO tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc > /dev/null
    if [ ! -f /etc/apt/sources.list.d/cran.list ]; then
        echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${VERSION_CODENAME}-cran40/" \
            | $SUDO tee /etc/apt/sources.list.d/cran.list > /dev/null
    fi
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

## GCC on Ubuntu. GCC_PREFER (ci/toolchain.env) is a PIN, not a wish: if it
## is set, install exactly that version or stop the build. Silently landing
## on whatever the PPA happened to have is how a bundle and the library
## inside it drift apart without anything saying so -- the baseline Linux
## bundle once carried a library built by GCC 14 inside a binary built by
## GCC 15, and the Windows DLL was built by MinGW 16.2.0 while inla.exe used
## 16.1.1, both discovered only after shipping. Mirrors the exact-pin-or-fail
## structure ci/deps-manylinux.sh already uses for GCC_TOOLSET.
##
## The version installed is EXPORTED (CC/CXX/FC), not just echoed, so
## ci/fetch-stiles.sh's BUILDINFO check and ci/build.sh's actual compile both
## use this exact compiler instead of each re-discovering "the newest gcc-N
## on disk" independently -- two independent discoveries can disagree even
## when both succeed. Both already prefer an inherited CXX (only falling
## back to their own probe if it is unset), so exporting it here is enough;
## nothing downstream needs to change.
if [ "${ID:-}" = "ubuntu" ]; then
    $SUDO apt-get install -y --no-install-recommends software-properties-common
    $SUDO add-apt-repository -y ppa:ubuntu-toolchain-r/test 2>/dev/null || true
    $SUDO apt-get update || true
    _R=$(cd "$(dirname "$0")/.." && pwd); [ -f "$_R/ci/toolchain.env" ] && . "$_R/ci/toolchain.env"

    _install_gcc() {
        $SUDO apt-get install -y --no-install-recommends \
            "gcc-$1" "g++-$1" "gfortran-$1" 2>/dev/null
    }
    _export_gcc() {
        echo "toolchain: gcc-$1"
        if [ -n "${GITHUB_ENV:-}" ]; then
            { echo "CC=gcc-$1"; echo "CXX=g++-$1"; echo "FC=gfortran-$1"; } >> "$GITHUB_ENV"
        fi
    }

    if [ -n "${GCC_PREFER:-}" ]; then
        if _install_gcc "$GCC_PREFER"; then
            _export_gcc "$GCC_PREFER"
        else
            echo "ERROR: gcc-$GCC_PREFER is not available (checked the" >&2
            echo "       ubuntu-toolchain-r/test PPA and the distro archive)." >&2
            echo "       Pinned by ci/toolchain.env; update it rather than falling" >&2
            echo "       back, so this project and sTiles keep using the same" >&2
            echo "       compiler. Set GCC_PREFER=\"\" there to allow a fallback." >&2
            exit 1
        fi
    else
        ## unpinned: take the newest available
        for v in 17 16 15 14; do
            if _install_gcc "$v"; then
                _export_gcc "$v"
                break
            fi
        done
    fi
fi
