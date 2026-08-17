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
ARMPL_VERSION=${ARMPL_VERSION:-26.07}
if ! ls -d /opt/arm/armpl_* >/dev/null 2>&1; then
    URL=${ARMPL_URL:-https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_${ARMPL_VERSION}/arm-performance-libraries_${ARMPL_VERSION}_macOS.tgz}
    if curl -sfL -o /tmp/armpl.tgz "$URL"; then
        mkdir -p /tmp/armpl && tar xzf /tmp/armpl.tgz -C /tmp/armpl
        PKG=$(find /tmp/armpl -name '*.pkg' 2>/dev/null | head -1)
        DIR=$(find /tmp/armpl -maxdepth 2 -type d -name 'armpl_*' 2>/dev/null | head -1)
        if [ -n "$PKG" ]; then
            sudo installer -pkg "$PKG" -target /
        elif [ -n "$DIR" ]; then
            sudo mkdir -p /opt/arm && sudo cp -R "$DIR" /opt/arm/
        else
            echo "WARNING: unrecognised ARMPL archive layout:"; ls -R /tmp/armpl | head -20
        fi
    else
        echo "WARNING: ARMPL ${ARMPL_VERSION} download failed; the build will use R's Rblas/Rlapack"
    fi
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
fi
[ -f "$DEPS/lib/libRmath.a" ] || { echo "ERROR: standalone Rmath was not built"; exit 1; }

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
