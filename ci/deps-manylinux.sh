#!/usr/bin/env bash
## Build dependencies inside the manylinux_2_28 container (AlmaLinux 8,
## glibc 2.28) used by the portable-binary CI job. Not meant for ordinary
## machines: use ci/deps-ubuntu.sh or ci/deps-fedora.sh there.
set -e

## Shared versions come from one file the maintainer owns, so this project and
## sTiles build with the same compiler rather than each picking its own.
_here=$(cd "$(dirname "$0")" && pwd)
[ -r "$_here/toolchain.env" ] && . "$_here/toolchain.env"

dnf -y install dnf-plugins-core epel-release
dnf config-manager --set-enabled powertools 2>/dev/null \
    || dnf config-manager --set-enabled crb 2>/dev/null || true

dnf -y install \
    gcc gcc-c++ gcc-gfortran make git-core findutils diffutils rsync patchelf \
    openblas-devel \
    gsl-devel \
    metis-devel \
    openssl-devel \
    suitesparse-devel \
    zlib-devel \
    libRmath-devel \
    R-core-devel \
    numactl-devel hwloc-devel libtool-ltdl-devel

## R-core-devel above provides only the HEADERS the external model packages
## compile against; the built binary uses dlopen mode, so the shipped
## bundle still needs no R anywhere.

## Eigen: the external packages use Eigen 3.4 API and EL8 packages 3.3 --
## vendor current 3.4 headers instead (header-only).
if [ ! -d /usr/local/include/eigen3/Eigen ]; then
    git clone --depth 1 --branch 3.4 https://gitlab.com/libeigen/eigen.git /tmp/eigen
    mkdir -p /usr/local/include/eigen3
    cp -r /tmp/eigen/Eigen /tmp/eigen/unsupported /usr/local/include/eigen3/
fi

## muParser: the package name varies across repos and EPEL8 may not carry
## it at all -- try both names, else build the (small, cmake-based) library
## from source into the container's system paths.
if ! dnf -y install muParser-devel 2>/dev/null \
   && ! dnf -y install muparser-devel 2>/dev/null; then
    dnf -y install cmake
    git clone --depth 1 https://github.com/beltoforion/muparser /tmp/muparser
    cmake -S /tmp/muparser -B /tmp/muparser/build \
          -DCMAKE_BUILD_TYPE=Release -DENABLE_SAMPLES=OFF -DENABLE_OPENMP=OFF \
          -DCMAKE_INSTALL_PREFIX=/usr
    cmake --build /tmp/muparser/build -j"$(nproc)"
    cmake --install /tmp/muparser/build
    ldconfig
fi

## Newest gcc-toolset available. These are recent GCCs that target the
## container's own glibc, so the binary keeps its low glibc floor while
## being compiled by a modern compiler. The image ships one already; take
## a newer one when the repositories offer it.
## GCC_TOOLSET (ci/toolchain.env) names the ONE toolset both this project and
## sTiles use, so a bundle and the library inside it are built by the same
## compiler. Taking "the newest available" independently on each side is how
## they drifted apart: this lane installed 15 while sTiles used the image's 14.
_want=${GCC_TOOLSET:-}
if [ -n "$_want" ]; then
    if dnf -y install "gcc-toolset-$_want" "gcc-toolset-$_want-gcc-c++" \
                      "gcc-toolset-$_want-gcc-gfortran" 2>/dev/null; then
        echo "toolchain: gcc-toolset-$_want (pinned by GCC_TOOLSET)"
    else
        echo "ERROR: gcc-toolset-$_want is not available in this container." >&2
        echo "       Pinned by ci/toolchain.env; update it rather than falling back," >&2
        echo "       so this project and sTiles keep using the same compiler." >&2
        exit 1
    fi
else
    for v in 17 16 15 14; do
        if dnf -y install "gcc-toolset-$v" "gcc-toolset-$v-gcc-c++" \
                          "gcc-toolset-$v-gcc-gfortran" 2>/dev/null; then
            echo "toolchain: gcc-toolset-$v (unpinned)"
            break
        fi
    done
fi

## Static archives for the libraries that get embedded into the binary.
dnf -y install openblas-static 2>/dev/null || true

## SIMDE is header-only and has no EL8 package: vendor the headers.
if [ ! -d /usr/local/include/simde ]; then
    git clone --depth 1 https://github.com/simd-everywhere/simde /tmp/simde-src
    cp -r /tmp/simde-src/simde /usr/local/include/simde
fi
