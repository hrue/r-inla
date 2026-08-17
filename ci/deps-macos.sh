#!/usr/bin/env bash
## Build environment for inla on macOS (Apple Silicon): Homebrew toolchain
## and libraries, plus the official CRAN R (the framework layout users
## have, including libR/libRmath/libRblas/libRlapack). metis has no
## reliable Homebrew formula, so it is built static from source.
set -e

## gnu-sed: gmrflib's dependency-file rule uses GNU sed syntax, which BSD
## sed silently mangles into an unparsable .d file.
brew install gcc gsl muparser eigen simde libtool openssl@3 cmake gnu-sed 2>/dev/null || true

## Official CRAN R for Apple Silicon (the user-realistic install).
if [ ! -d /Library/Frameworks/R.framework ]; then
    PKG=$(curl -s https://cran.r-project.org/bin/macosx/big-sur-arm64/base/ \
          | grep -oE 'R-[0-9.]+-arm64\.pkg' | sort -V | tail -1)
    curl -s -o /tmp/R.pkg "https://cran.r-project.org/bin/macosx/big-sur-arm64/base/$PKG"
    sudo installer -pkg /tmp/R.pkg -target /
fi

## Arm Performance Libraries: the BLAS/LAPACK the upstream macOS binaries
## use. Installed from Arm's public download; if the URL has moved (version
## bumps), the build falls back to R's own Rblas/Rlapack automatically.
ARMPL_VERSION=${ARMPL_VERSION:-24.10}
if ! ls -d /opt/arm/armpl_* >/dev/null 2>&1; then
    URL="https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_${ARMPL_VERSION}/arm-performance-libraries_${ARMPL_VERSION}_macOS.dmg"
    if curl -sfL -o /tmp/armpl.dmg "$URL"; then
        hdiutil attach -quiet -nobrowse -mountpoint /Volumes/armpl /tmp/armpl.dmg
        PKG=$(ls /Volumes/armpl/*.pkg 2>/dev/null | head -1)
        [ -n "$PKG" ] && sudo installer -pkg "$PKG" -target /
        hdiutil detach -quiet /Volumes/armpl || true
    else
        echo "WARNING: ARMPL ${ARMPL_VERSION} download failed; the build will use R's Rblas/Rlapack"
    fi
fi

DEPS=$HOME/inla-macos-deps
mkdir -p "$DEPS/lib" "$DEPS/include"

if [ ! -f "$DEPS/lib/libmetis.a" ]; then
    ## The scivision mirror carries a modern cmake build that also handles
    ## GKlib; the policy flag keeps any old cmake_minimum_required lines
    ## acceptable to current CMake.
    git clone --depth 1 https://github.com/scivision/METIS /tmp/metis
    cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/metis/build -j"$(sysctl -n hw.ncpu)"
    cmake --install /tmp/metis/build
fi

echo "OK: macOS deps ready (brew + CRAN R + $DEPS)"
