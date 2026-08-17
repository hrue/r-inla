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
    eigen3-devel \
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
      ## the RNG stub's filename has changed across R versions: take
      ## whatever .c the standalone dir provides (minus its test program)
      STD=$(ls standalone/*.c 2>/dev/null | grep -v test || true)
      for f in *.c $STD; do
          $MINGW_CC -O2 -DMATHLIB_STANDALONE \
              -I. -I"$DEPS/R-win/include" -c "$f" -o "$(basename "$f" .c).o"
      done
      $TRIPLET-ar rcs "$DEPS/lib/libRmath.a" ./*.o )
fi

## ---- dlfcn-win32: the Windows port of dlfcn.h, which the external model
## packages include. Fedora packages it as mingw64-dlfcn/ucrt64-dlfcn; build
## it from source when neither is available.
if ! $SUDO dnf -y install ucrt64-dlfcn 2>/dev/null \
   && ! $SUDO dnf -y install mingw64-dlfcn 2>/dev/null \
   && [ ! -f "$DEPS/lib/libdl.a" ]; then
    git clone --depth 1 https://github.com/dlfcn-win32/dlfcn-win32 /tmp/dlfcn
    cmake -S /tmp/dlfcn -B /tmp/dlfcn/build \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/dlfcn/build -j"$(nproc)"
    cmake --install /tmp/dlfcn/build
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
    ## The scivision mirror carries a modern cmake build that also handles
    ## GKlib and cross-compilation.
    git clone --depth 1 https://github.com/scivision/METIS /tmp/metis
    cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX="$DEPS"
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

## ---- libltdl: the UCRT sysroot has no ltdl.h, so take the header and the
## import library from the mingw64 sysroot into our own deps tree (upstream
## links that same library into UCRT builds). --------------------------------
if [ ! -f "$DEPS/include/ltdl.h" ]; then
    H=$(find /usr/x86_64-w64-mingw32*/sys-root/mingw/include -name ltdl.h 2>/dev/null | head -1 || true)
    [ -n "$H" ] || { echo "ERROR: ltdl.h not found in any mingw sysroot"; exit 1; }
    cp "$H" "$DEPS/include/"
    D=$(dirname "$H")/libltdl
    [ -d "$D" ] && cp -r "$D" "$DEPS/include/"
    L=$(find /usr/x86_64-w64-mingw32*/sys-root/mingw/lib -name 'libltdl.dll.a' 2>/dev/null | head -1 || true)
    [ -n "$L" ] && cp "$L" "$DEPS/lib/"
fi

## ---- SIMDE headers (header-only, host-independent) --------------------------
## From upstream rather than the distribution package: Fedora's simde
## expects hedley.h from a separate package, while the upstream tree keeps
## it inside simde/, so the copied headers stay self-contained.
if [ ! -d "$DEPS/include/simde" ]; then
    git clone --depth 1 https://github.com/simd-everywhere/simde /tmp/simde-src
    cp -r /tmp/simde-src/simde "$DEPS/include/simde"
fi

## ---- Eigen headers (header-only, so the host package serves the cross
## build as well) -------------------------------------------------------------
if [ ! -d "$DEPS/include/eigen3" ]; then
    for d in /usr/include/eigen3 /usr/x86_64-w64-mingw32ucrt/sys-root/mingw/include/eigen3 \
             /usr/x86_64-w64-mingw32/sys-root/mingw/include/eigen3; do
        [ -d "$d" ] && { cp -r "$d" "$DEPS/include/eigen3"; break; }
    done
fi
[ -d "$DEPS/include/eigen3/Eigen" ] || { echo "ERROR: no Eigen headers found"; exit 1; }

echo "OK: mingw deps ready under $DEPS"
