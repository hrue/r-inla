#!/usr/bin/env bash
## Cross-compilation environment for the Windows inla.exe, inside a Fedora
## container: the ucrt64-* MinGW toolchain plus every library the link
## needs, built static where no UCRT package exists. Also unpacks the
## official Windows R distribution (for its headers) and cross-builds the
## standalone Rmath library from the R sources (the code uses
## MATHLIB_STANDALONE, whose unprefixed symbols R.dll does not export).
set -e

## Versions are pinned in one file the maintainer owns.
_ROOT=$(cd "$(dirname "$0")/.." && pwd)
[ -f "$_ROOT/ci/toolchain.env" ] && . "$_ROOT/ci/toolchain.env"

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
## Prefer a UCRT libltdl when one is packaged; the mingw64 build installed
## above is the fallback (its DLL works, being a thin dlopen wrapper).
$SUDO dnf -y install ucrt64-libltdl 2>/dev/null || true
$SUDO dnf -y install ucrt64-openssl 2>/dev/null || $SUDO dnf -y install mingw64-openssl 2>/dev/null || true

## OpenMP runtime for the cross compiler. INLA requires OpenMP (the sources
## refuse to compile without it), and the link needs libgomp.spec next to
## the target's libraries. Package naming varies, so try the candidates and
## report the toolchain's actual contents when none of them provides it.
for p in ucrt64-libgomp mingw64-libgomp ucrt64-gcc-libgomp ucrt64-openmp; do
    $SUDO dnf -y install "$p" 2>/dev/null && break
done
if ! find /usr/lib/gcc/$TRIPLET /usr/$TRIPLET/sys-root/mingw/lib \
          -name 'libgomp.spec' 2>/dev/null | grep -q .; then
    echo "WARNING: no libgomp.spec for $TRIPLET -- OpenMP link will fail"
    echo "  gcc lib dir:"; ls -R /usr/lib/gcc/$TRIPLET 2>/dev/null | head -30
    echo "  candidate packages:"
    dnf -q list --available '*gomp*' 2>/dev/null | grep -iE 'ucrt|mingw' | head -10 || true
    echo "  files providing libgomp.spec:"
    dnf -q provides '*/libgomp.spec' 2>/dev/null | head -10 || true
fi

DEPS=/opt/mingw-deps
$SUDO mkdir -p "$DEPS"/lib "$DEPS"/include
$SUDO chmod 777 "$DEPS"

## ---- Windows R: HEADERS ONLY ----------------------------------------------
## Nothing links or ships a DLL out of this tree. It is unpacked for
## include/ (Rmath.h0 and friends) and nothing else.
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

    ## The matching header, generated from R's own template as the
    ## standalone build does (the shipped header is written for use from
    ## inside R and omits part of the standalone API).
    RVER=$(echo "$R_SRC" | sed -e 's/^R-//' -e 's/\.tar\.gz$//')
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
    grep -q 'double.*lbeta' "$DEPS/include/Rmath.h" \
        || { echo "ERROR: generated Rmath.h lacks lbeta"; exit 1; }
fi

## ---- OpenBLAS: this build's BLAS and LAPACK --------------------------------
## Built from source with THIS toolchain and linked statically, so no BLAS
## library reaches the shipped bundle.
##
## It replaces R's own Rblas.dll/Rlapack.dll, which cannot be used here. The
## exe imports those by name, they import R.dll, so the bundle then has to
## ship R.dll beside the exe; Windows loads that copy at process start, while
## rgeneric later dlopens the user's R.dll through R_HOME. A process holding
## two R.dll copies segfaults the moment the embedded R starts, even when
## both copies are the same version. Linux and macOS never reach that state
## because their BLAS is MKL or Accelerate, so no R library is ever
## preloaded, which is why rgeneric worked everywhere except here.
##
## Single threaded on purpose: INLA parallelises above the BLAS with OpenMP,
## and the reference Rblas this replaces was single threaded too, so the
## threading behaviour does not change while the kernels get much faster.
##
## Runtime dispatch, not a fixed TARGET. This is what makes the Windows
## binary both fast and portable, and the two are not in tension here:
## OpenBLAS picks its kernels from cpuid at startup, so one build runs the
## SKYLAKEX or ZEN path on a modern machine and still starts on an old one.
##
## TARGET sets the BASELINE that a machine matching nothing else falls back
## to. PRESCOTT is plain x86-64, which is exactly the floor the rest of
## inla.exe has (it compiles with -mtune=generic and no -march), so this adds
## no CPU requirement at all. An earlier version pinned TARGET=NEHALEM with no
## dispatch, which was both slower on modern hardware AND raised the floor to
## roughly 2009 machines: the worst of both.
##
## DYNAMIC_LIST keeps the cost bounded. Plain DYNAMIC_ARCH builds every kernel
## OpenBLAS knows, on every run, since nothing caches $DEPS between runs; this
## list covers the hardware anyone actually runs INLA on. Add a name here
## rather than dropping back to a fixed TARGET if a gap shows up.
if [ ! -f "$DEPS/lib/libopenblas.a" ]; then
    OB=${OPENBLAS_VERSION:-0.3.29}
    rm -rf /tmp/openblas
    git clone -q --depth 1 --branch "v$OB" \
        https://github.com/OpenMathLib/OpenBLAS /tmp/openblas
    ## CROSS=1 stops the build from running the test binaries it produces,
    ## which are Windows executables here.
    if ! make -C /tmp/openblas -j"$(nproc)" \
            HOSTCC=gcc CC=$MINGW_CC FC=$TRIPLET-gfortran \
            AR=$TRIPLET-ar RANLIB=$TRIPLET-ranlib \
            CROSS=1 BINARY=64 \
            DYNAMIC_ARCH=1 TARGET=PRESCOTT \
            DYNAMIC_LIST="NEHALEM SANDYBRIDGE HASWELL SKYLAKEX ZEN" \
            NO_SHARED=1 USE_THREAD=0 NO_LAPACKE=1 \
            > /tmp/openblas.log 2>&1; then
        echo "ERROR: OpenBLAS $OB failed to cross-build"
        tail -40 /tmp/openblas.log
        exit 1
    fi
    ## the same vars as the build: install derives the library name from
    ## them, and a mismatch makes it look for an archive that is not there
    make -C /tmp/openblas PREFIX="$DEPS" NO_SHARED=1 \
         DYNAMIC_ARCH=1 TARGET=PRESCOTT install >/dev/null
    ## some versions install the archive under its target-specific name only
    if [ ! -f "$DEPS/lib/libopenblas.a" ]; then
        A=$(ls "$DEPS"/lib/libopenblas*.a 2>/dev/null | head -1 || true)
        [ -n "$A" ] && cp -f "$A" "$DEPS/lib/libopenblas.a"
    fi
fi
[ -f "$DEPS/lib/libopenblas.a" ] \
    || { echo "ERROR: no libopenblas.a under $DEPS/lib"; ls "$DEPS/lib"; exit 1; }

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

## ---- mimalloc: fast allocator, as in the upstream Windows recipe -----------
## Upstream links libmimalloc.dll.a first in EXTLIBS2 plus
## -Wl,--undefined=mi_version, which routes the process allocations through
## mimalloc. GMRFLib is allocation-heavy, so this is a real speedup on
## Windows, where the UCRT allocator is comparatively slow.
if [ ! -f "$DEPS/lib/libmimalloc.dll.a" ]; then
    git clone --depth 1 --branch v2.2.4 https://github.com/microsoft/mimalloc /tmp/mimalloc \
      || git clone --depth 1 https://github.com/microsoft/mimalloc /tmp/mimalloc
    cmake -S /tmp/mimalloc -B /tmp/mimalloc/build \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC \
        -DCMAKE_CXX_COMPILER=${MINGW_CC%-gcc}-g++ \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON \
        -DMI_BUILD_STATIC=OFF -DMI_BUILD_TESTS=OFF -DMI_BUILD_OBJECT=OFF \
        -DMI_WIN_REDIRECT=OFF \
        -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/mimalloc/build -j"$(nproc)"
    cmake --install /tmp/mimalloc/build
    ## normalize install layout: some versions install under lib/mimalloc-2.x
    find "$DEPS" -name 'libmimalloc*.dll.a' -exec cp -n {} "$DEPS/lib/libmimalloc.dll.a" \; 2>/dev/null || true
    find "$DEPS" -name 'libmimalloc*.dll' -exec cp -n {} "$DEPS/lib/" \; 2>/dev/null || true
    ls "$DEPS/lib"/libmimalloc* 
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
    ## GitHub throttles anonymous downloads (cmake's FetchContent pulls
    ## GKlib from there), so retry rather than fail the whole lane.
    for attempt in 1 2 3; do
        cmake -S /tmp/metis -B /tmp/metis/build \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
            -DCMAKE_INSTALL_PREFIX="$DEPS" && break
        echo "METIS configure failed (attempt $attempt/3), retrying"
        sleep $((attempt * 20))
    done
    cmake --build /tmp/metis/build -j"$(nproc)"
    cmake --install /tmp/metis/build
fi

## muParser as a DLL, as the upstream Windows build links it: inla's eval.c
## uses muParser's C API, whose header declares the entry points
## __declspec(dllimport) on Windows, so the symbols come out as __imp_mup*
## and only an import library resolves them.
if [ ! -f "$DEPS/lib/libmuparser.dll.a" ]; then
    git clone --depth 1 https://github.com/beltoforion/muparser /tmp/muparser
    cmake -S /tmp/muparser -B /tmp/muparser/build \
        -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=$MINGW_CC -DCMAKE_CXX_COMPILER=$MINGW_CXX \
        -DCMAKE_BUILD_TYPE=Release -DENABLE_SAMPLES=OFF -DENABLE_OPENMP=OFF \
        -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$DEPS"
    cmake --build /tmp/muparser/build -j"$(nproc)"
    cmake --install /tmp/muparser/build
    ## cmake puts the DLL in bin/ and the import library in lib/
    cp -f "$DEPS"/bin/*muparser*.dll "$DEPS/lib/" 2>/dev/null || true
fi
ls "$DEPS"/lib/libmuparser.dll.a >/dev/null 2>&1 \
    || { echo "ERROR: no muParser import library was produced"; ls "$DEPS"/lib "$DEPS"/bin; exit 1; }

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
