#!/usr/bin/env bash
## Cross-compilation environment for the Windows inla.exe, inside a Fedora
## container: the ucrt64-* MinGW toolchain plus every library the link
## needs, built static where no UCRT package exists. Also unpacks the
## official Windows R distribution (for R.dll, Rblas/Rlapack, headers) and
## cross-builds the standalone Rmath library from the R sources (the code
## uses MATHLIB_STANDALONE, whose unprefixed symbols R.dll does not export).
set -e

TRIPLET=x86_64-w64-mingw32ucrt
MINGW_CC=$TRIPLET-gcc
MINGW_CXX=$TRIPLET-g++

SUDO=""
[ "$(id -u)" != 0 ] && SUDO=sudo

$SUDO dnf -y install \
    ucrt64-gcc ucrt64-gcc-c++ ucrt64-gcc-gfortran ucrt64-winpthreads \
    mingw64-eigen3 mingw64-libltdl \
    simde-devel \
    git-core make cmake findutils diffutils rsync wget innoextract zip \
    gcc gcc-c++ R-core-devel libRmath-devel

## Optional UCRT packages: use them when present, otherwise the source
## builds below cover the gaps.
$SUDO dnf -y install ucrt64-zlib 2>/dev/null || true
$SUDO dnf -y install ucrt64-openssl 2>/dev/null || $SUDO dnf -y install mingw64-openssl 2>/dev/null || true

DEPS=/opt/mingw-deps
$SUDO mkdir -p "$DEPS"/lib "$DEPS"/include
$SUDO chmod 777 "$DEPS"

## ---- Windows R (R.dll, Rblas.dll, Rlapack.dll, include/) ------------------
if [ ! -f "$DEPS/R-win/bin/x64/R.dll" ]; then
    R_EXE=$(wget -qO- https://cran.r-project.org/bin/windows/base/ \
            | grep -oE 'R-[0-9.]+-win\.exe' | head -1)
    wget -q "https://cran.r-project.org/bin/windows/base/$R_EXE" -O /tmp/R-win.exe
    innoextract -s -d /tmp/R-win /tmp/R-win.exe
    mv /tmp/R-win/app "$DEPS/R-win"
fi

## ---- standalone Rmath, cross-compiled static -------------------------------
if [ ! -f "$DEPS/lib/libRmath.a" ]; then
    R_SRC=$(wget -qO- https://cran.r-project.org/src/base/R-4/ \
            | grep -oE 'R-[0-9.]+\.tar\.gz' | sort -V | tail -1)
    wget -q "https://cran.r-project.org/src/base/R-4/$R_SRC" -O /tmp/R-src.tar.gz
    mkdir -p /tmp/R-src && tar xzf /tmp/R-src.tar.gz -C /tmp/R-src --strip-components=1
    ( cd /tmp/R-src/src/nmath
      ## minimal configuration header the sources expect from configure
      cat > config.h <<'EOF'
#define HAVE_EXPM1 1
#define HAVE_HYPOT 1
#define HAVE_LOG1P 1
#define HAVE_WORKING_LOG1P 1
EOF
      for f in *.c standalone/std_unif.c; do
          $MINGW_CC -O2 -DMATHLIB_STANDALONE \
              -I. -I"$DEPS/R-win/include" -c "$f" -o "$(basename "$f" .c).o"
      done
      $TRIPLET-ar rcs "$DEPS/lib/libRmath.a" ./*.o )
fi

## ---- gsl, metis, muparser: static mingw builds ------------------------------
if [ ! -f "$DEPS/lib/libgsl.a" ]; then
    GSL=$(wget -qO- https://ftp.gnu.org/gnu/gsl/ | grep -oE 'gsl-2\.[0-9.]+\.tar\.gz' | sort -V | tail -1)
    wget -q "https://ftp.gnu.org/gnu/gsl/$GSL" -O /tmp/gsl.tar.gz
    mkdir -p /tmp/gsl && tar xzf /tmp/gsl.tar.gz -C /tmp/gsl --strip-components=1
    ( cd /tmp/gsl
      ./configure --host=$TRIPLET --prefix="$DEPS" --disable-shared --enable-static --quiet
      make -j"$(nproc)" >/dev/null && make install >/dev/null )
fi

if [ ! -f "$DEPS/lib/libmetis.a" ]; then
    git clone --depth 1 https://github.com/KarypisLab/METIS /tmp/metis 2>/dev/null || \
        wget -q http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz -O /tmp/metis.tar.gz
    if [ ! -d /tmp/metis ]; then mkdir -p /tmp/metis && tar xzf /tmp/metis.tar.gz -C /tmp/metis --strip-components=1; fi
    cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC \
        -DCMAKE_BUILD_TYPE=Release -DSHARED=OFF -DCMAKE_INSTALL_PREFIX="$DEPS" \
        -DGKLIB_PATH=/tmp/metis/GKlib 2>/dev/null || true
    cmake --build /tmp/metis/build -j"$(nproc)"
    cmake --install /tmp/metis/build
fi

if [ ! -f "$DEPS/lib/libmuparser.a" ]; then
    git clone --depth 1 https://github.com/beltoforion/muparser /tmp/muparser
    cmake -S /tmp/muparser -B /tmp/muparser/build \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC -DCMAKE_CXX_COMPILER=$MINGW_CXX \
        -DCMAKE_BUILD_TYPE=Release -DENABLE_SAMPLES=OFF -DENABLE_OPENMP=OFF \
        -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/muparser/build -j"$(nproc)"
    cmake --install /tmp/muparser/build
fi

## ---- SIMDE headers in a clean include dir (they are host-independent) ------
if [ ! -d "$DEPS/include/simde" ]; then
    cp -r /usr/include/simde "$DEPS/include/simde"
fi

echo "OK: mingw deps ready under $DEPS"
