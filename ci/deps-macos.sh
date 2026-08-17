#!/usr/bin/env bash
## Build environment for inla on macOS (Apple Silicon): Homebrew toolchain
## and libraries, plus the official CRAN R (the framework layout users
## have, including libR/libRmath/libRblas/libRlapack). metis has no
## reliable Homebrew formula, so it is built static from source.
set -e

brew install gcc gsl muparser eigen simde libtool openssl@3 cmake 2>/dev/null || true

## Official CRAN R for Apple Silicon (the user-realistic install).
if [ ! -d /Library/Frameworks/R.framework ]; then
    PKG=$(curl -s https://cran.r-project.org/bin/macosx/big-sur-arm64/base/ \
          | grep -oE 'R-[0-9.]+-arm64\.pkg' | sort -V | tail -1)
    curl -s -o /tmp/R.pkg "https://cran.r-project.org/bin/macosx/big-sur-arm64/base/$PKG"
    sudo installer -pkg /tmp/R.pkg -target /
fi

DEPS=$HOME/inla-macos-deps
mkdir -p "$DEPS/lib" "$DEPS/include"

if [ ! -f "$DEPS/lib/libmetis.a" ]; then
    curl -sL -o /tmp/metis.tar.gz \
        https://github.com/KarypisLab/METIS/archive/refs/tags/v5.2.1.tar.gz 2>/dev/null || true
    if [ -s /tmp/metis.tar.gz ]; then
        mkdir -p /tmp/metis && tar xzf /tmp/metis.tar.gz -C /tmp/metis --strip-components=1
    else
        git clone --depth 1 https://github.com/scivision/METIS /tmp/metis
    fi
    ## METIS >= 5.2 needs GKlib separately; the scivision mirror bundles a
    ## cmake build that handles both.
    cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/metis/build -j"$(sysctl -n hw.ncpu)"
    cmake --install /tmp/metis/build
fi

echo "OK: macOS deps ready (brew + CRAN R + $DEPS)"
