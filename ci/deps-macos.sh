#!/usr/bin/env bash
## Build environment for inla on macOS (Apple Silicon): Homebrew toolchain
## and libraries, plus the official CRAN R (the framework layout users
## have, including libR/libRmath/libRblas/libRlapack). metis has no
## reliable Homebrew formula, so it is built static from source.
set -e

## gnu-sed: gmrflib's dependency-file rule uses GNU sed syntax, which BSD
## sed silently mangles into an unparsable .d file.
## dylibbundler makes the built binary self-contained (ci/package-macos.sh).
brew install gcc gsl muparser eigen simde libtool openssl@3 cmake gnu-sed dylibbundler 2>/dev/null || true

## Official CRAN R for Apple Silicon (the user-realistic install).
if [ ! -d /Library/Frameworks/R.framework ]; then
    if [ "$(uname -m)" = arm64 ]; then
        BASE=https://cran.r-project.org/bin/macosx/big-sur-arm64/base
        PAT='R-[0-9.]+-arm64\.pkg'
    else
        BASE=https://cran.r-project.org/bin/macosx/big-sur-x86_64/base
        PAT='R-[0-9.]+-x86_64\.pkg'
    fi
    PKG=$(curl -s "$BASE/" | grep -oE "$PAT" | sort -V | tail -1)
    [ -n "$PKG" ] || { echo "ERROR: no R installer found at $BASE"; exit 1; }
    curl -s -o /tmp/R.pkg "$BASE/$PKG"
    sudo installer -pkg /tmp/R.pkg -target /
fi

## Arm Performance Libraries: the BLAS/LAPACK the upstream macOS binaries
## use. Installed from Arm's public download; if the URL has moved (version
## bumps), the build falls back to R's own Rblas/Rlapack automatically.
## Apple Silicon only: the Intel lane uses the serial OpenBLAS built below.
## Homebrew's cask installs ARMPL under /opt/arm, which is simpler and more
## reliable than fetching from Arm's CDN. Retried: that CDN drops
## connections often enough to take out a whole run.
if [ "$(uname -m)" = arm64 ] && ! ls -d /opt/arm/armpl_* >/dev/null 2>&1; then
    for attempt in 1 2 3; do
        brew install --cask arm-performance-libraries && break
        echo "ARMPL download failed (attempt $attempt/3), retrying"
        sleep $((attempt * 20))
    done
    ls -d /opt/arm/armpl_* >/dev/null 2>&1 \
        || echo "WARNING: ARMPL unavailable; the build will use R's Rblas/Rlapack"
fi

DEPS=$HOME/inla-macos-deps
mkdir -p "$DEPS/lib" "$DEPS/include"

## Standalone Rmath: the CRAN framework ships libR/libRblas/libRlapack but
## not the standalone math library the sources use (rmath.h defines
## MATHLIB_STANDALONE), so build it from R's own sources.
if [ ! -f "$DEPS/lib/libRmath.a" ]; then
    RHOME=$(R RHOME)
    RVER=$("$RHOME/bin/R" --version | head -1 | awk '{print $3}')
    curl -sL -o /tmp/R-src.tar.gz "https://cran.r-project.org/src/base/R-4/R-${RVER}.tar.gz"
    mkdir -p /tmp/R-src && tar xzf /tmp/R-src.tar.gz -C /tmp/R-src --strip-components=1
    ( cd /tmp/R-src/src/nmath
      cat > config.h <<'EOF'
#define HAVE_EXPM1 1
#define HAVE_HYPOT 1
#define HAVE_LOG1P 1
#define HAVE_WORKING_LOG1P 1
EOF
      CCB=$(ls "$(brew --prefix)"/bin/gcc-1[0-9] 2>/dev/null | sort -V | tail -1 || echo cc)
      STD=$(ls standalone/*.c 2>/dev/null | grep -v test || true)
      for f in *.c $STD; do
          "$CCB" -O2 -DMATHLIB_STANDALONE -I. -I"$RHOME/include" \
                 -c "$f" -o "$(basename "$f" .c).o"
      done
      ar rcs "$DEPS/lib/libRmath.a" ./*.o )

    ## The matching header, generated from R's own template the way the
    ## standalone build does. The installed framework header is written for
    ## use from inside R and does not declare everything the standalone
    ## library provides (lbeta, the M_LN_SQRT_* constants).
    ## configure generates Rmath.h0 from Rmath.h0.in, and a source tarball
    ## carries only the template, so accept either.
    TPL=$(ls /tmp/R-src/src/include/Rmath.h0 /tmp/R-src/src/include/Rmath.h0.in 2>/dev/null | head -1 || true)
    [ -n "$TPL" ] || { echo "ERROR: no Rmath.h0 template in the R sources"; exit 1; }
    ## configure normally fills these in; do the same substitutions here.
    ## HAVE_WORKING_LOG1P matters: without it the header falls into a block
    ## that declares Rlog1p inside extern "C", which is not valid C.
    sed -e "s/@PACKAGE_VERSION@/$RVER/g" \
        -e 's|@RMATH_HAVE_WORKING_LOG1P@|#define HAVE_WORKING_LOG1P 1|' \
        -e 's|^#undef MATHLIB_STANDALONE|#define MATHLIB_STANDALONE 1|' \
        -e 's|^/\* #undef MATHLIB_STANDALONE \*/|#define MATHLIB_STANDALONE 1|' \
        -e 's/@[A-Za-z0-9_]*@//g' \
        "$TPL" > "$DEPS/include/Rmath.h"
    if grep -q '@[A-Za-z0-9_]*@' "$DEPS/include/Rmath.h"; then
        echo "ERROR: unsubstituted tokens remain in the generated Rmath.h"
        grep -n '@[A-Za-z0-9_]*@' "$DEPS/include/Rmath.h" | head
        exit 1
    fi
fi
[ -f "$DEPS/lib/libRmath.a" ] || { echo "ERROR: standalone Rmath was not built"; exit 1; }
[ -f "$DEPS/include/Rmath.h" ] || { echo "ERROR: standalone Rmath.h was not generated"; exit 1; }
grep -q 'double.*lbeta' "$DEPS/include/Rmath.h" || { echo "ERROR: generated Rmath.h lacks lbeta"; exit 1; }

if [ ! -f "$DEPS/lib/libmetis.a" ]; then
    ## The scivision mirror carries a modern cmake build that also handles
    ## GKlib; the policy flag keeps any old cmake_minimum_required lines
    ## acceptable to current CMake.
    git clone --depth 1 https://github.com/scivision/METIS /tmp/metis
    ## GitHub throttles anonymous downloads (cmake's FetchContent pulls
    ## GKlib from there), so retry rather than fail the whole lane.
    for attempt in 1 2 3; do
        cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
            -DCMAKE_INSTALL_PREFIX="$DEPS" && break
        echo "METIS configure failed (attempt $attempt/3), retrying"
        sleep $((attempt * 20))
    done
    cmake --build /tmp/metis/build -j"$(sysctl -n hw.ncpu)"
    cmake --install /tmp/metis/build
fi

echo "OK: macOS deps ready (brew + CRAN R + $DEPS)"

## Serial, locking-safe OpenBLAS for the Intel lane (arm64 uses ARMPL).
## Accelerate's threading cannot be controlled reliably, which breaks
## INLA's nested parallelism; this build never threads itself. Same
## recipe as the sTiles Intel lane, so both sides of the pairing sit on
## the same BLAS behavior.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
if [ "$(uname -m)" = x86_64 ] && [ ! -f "$ROOT/local/openblas/lib/libopenblas.dylib" ]; then
    for attempt in 1 2 3; do
        git clone --depth 1 --branch v0.3.29 \
            https://github.com/OpenMathLib/OpenBLAS.git /tmp/openblas-src && break
        echo "OpenBLAS clone failed (attempt $attempt/3), retrying"; sleep $((attempt*20)); rm -rf /tmp/openblas-src
    done
    make -C /tmp/openblas-src -j3 TARGET=HASWELL \
        USE_THREAD=0 USE_OPENMP=0 USE_LOCKING=1 \
        libs netlib shared > /tmp/ob.log 2>&1 || { tail -30 /tmp/ob.log; exit 1; }
    touch /tmp/openblas-src/lib.grd
    make -C /tmp/openblas-src install PREFIX="$ROOT/local/openblas" >> /tmp/ob.log 2>&1
    ls -l "$ROOT/local/openblas/lib/"libopenblas*.dylib
fi
