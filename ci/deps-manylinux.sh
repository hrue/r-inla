#!/usr/bin/env bash
## Build dependencies inside the manylinux_2_28 container (AlmaLinux 8,
## glibc 2.28) used by the portable-binary CI job. Not meant for ordinary
## machines: use ci/deps-ubuntu.sh or ci/deps-fedora.sh there.
set -e

dnf -y install dnf-plugins-core epel-release
dnf config-manager --set-enabled powertools 2>/dev/null \
    || dnf config-manager --set-enabled crb 2>/dev/null || true

dnf -y install \
    gcc gcc-c++ gcc-gfortran make git-core findutils diffutils rsync patchelf \
    openblas-devel \
    gsl-devel \
    metis-devel \
    muparser-devel \
    openssl-devel \
    suitesparse-devel \
    zlib-devel \
    eigen3-devel \
    libRmath-devel \
    numactl-devel hwloc-devel libtool-ltdl-devel

## SIMDE is header-only and has no EL8 package: vendor the headers.
if [ ! -d /usr/local/include/simde ]; then
    git clone --depth 1 https://github.com/simd-everywhere/simde /tmp/simde-src
    cp -r /tmp/simde-src/simde /usr/local/include/simde
fi
