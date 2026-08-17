#!/usr/bin/env bash
## Build inla on macOS (Apple Silicon) with the Homebrew GCC toolchain,
## following the same stages as ci/build.sh. Differences from Linux, per
## the upstream macOS recipe: ld64's -force_load instead of
## --whole-archive, no NUMA/hwloc/CLONE_TARGETS (Linux-only), BLAS/LAPACK
## from R's own libRblas/libRlapack, static libstdc++/libgcc, and rpaths
## into R's lib dir so the binary runs without environment setup.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
DEPS=${DEPS:-$HOME/inla-macos-deps}
JOBS=${JOBS:-$(sysctl -n hw.ncpu)}
BREW=$(brew --prefix)

## Homebrew's gcc: pick the newest versioned binaries.
CC=${CC:-$(ls "$BREW"/bin/gcc-1[0-9] 2>/dev/null | sort -V | tail -1 || true)}
CXX=${CXX:-$(ls "$BREW"/bin/g++-1[0-9] 2>/dev/null | sort -V | tail -1 || true)}
FC=${FC:-$(ls "$BREW"/bin/gfortran-1[0-9] 2>/dev/null | sort -V | tail -1 || true)}
[ -n "$CC" ] && [ -n "$CXX" ] && [ -n "$FC" ] || { echo "ERROR: Homebrew gcc not found"; exit 1; }

RHOME=$(R RHOME)
EPATH=$ROOT/external-packages
TAG=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo devel)
echo "== building $TAG for macOS arm64 with CC=$CC ($($CC --version | head -1)) =="

## Optimization: the upstream Apple Silicon configuration. -mcpu=apple-m1
## implies the armv8.5-a baseline upstream targets and adds the scheduling
## model for the same cores. LTO stays opt-in (LTO=1).
OPT=${OPT:-fast}
LTO=${LTO:-0}
case "$OPT" in
    fast) OPTFLAGS="-O3 -ftree-vectorize -funroll-loops -fvariable-expansion-in-unroller -fno-optimize-sibling-calls -ftracer" ;;
    safe) OPTFLAGS="-O2 -ftree-vectorize" ;;
    *)    echo "ERROR: OPT must be fast or safe"; exit 1 ;;
esac
[ "$LTO" = 1 ] && OPTFLAGS="$OPTFLAGS -flto=auto -ffat-lto-objects"
echo "== optimization: OPT=$OPT LTO=$LTO =="

## Oldest macOS the artifact must load on. Without this the deployment
## target defaults to the BUILD machine's OS, which would lock the binary
## out of every older Mac. 11.0 is the first macOS on Apple Silicon, so
## nothing is excluded on either architecture.
export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
echo "== macOS deployment target: $MACOSX_DEPLOYMENT_TARGET =="

## Apple Silicon gets the M-series tuning; Intel Macs the generic 64-bit
## baseline upstream uses (every Intel Mac on a current macOS is x86-64).
ARCH=$(uname -m)
if [ "$ARCH" = arm64 ]; then ARCHFLAGS="-mcpu=apple-m1"; else ARCHFLAGS="-mtune=generic -m64"; fi

FLAGS="$OPTFLAGS $ARCHFLAGS -pipe -pthread \
 -fopenmp -fopenmp-simd -flax-vector-conversions \
 -DINLA_WITH_SIMDE -DINLA_WITH_DEVEL -DINLA_WITH_CLONE_TARGETS \
 -DINLA_WITH_EXTERNAL_PACKAGES -DINLA_WITH_MUPARSER \
 -DGITCOMMIT=$TAG \
 -I$BREW/include -I$DEPS/include"

mkdir -p "$PREFIX"/bin "$PREFIX"/lib "$PREFIX"/include "$PREFIX/include.boot"
ln -sfn "$ROOT/gmrflib" "$PREFIX/include.boot/GMRFLib"
FLAGS="$FLAGS -I$PREFIX/include.boot"

## R linkage, as on Linux: 1 links libR at build time, 2 loads the running
## machine's libR through libltdl when an rgeneric model first appears (so
## the binary starts with any R or none). -L$DEPS/lib supplies the
## standalone Rmath built by ci/deps-macos.sh: the framework ships libR but
## not that one.
WITH_LIBR=${WITH_LIBR:-1}
if [ "$WITH_LIBR" = 2 ]; then
    ## The R include path stays: libR is not linked, but rmath.h still
    ## includes <Rmath.h> for the standalone math library.
    RLIB_INC="-DINLA_WITH_LIBR -DINLA_WITH_LIBR_DLOPEN -I$RHOME/include"
    RLIB_LIB="-L$DEPS/lib -lRmath"
else
    RLIB_INC="-DINLA_WITH_LIBR -I$RHOME/include"
    RLIB_LIB="-L$RHOME/lib -L$DEPS/lib -lR -lRmath"
fi
echo "== R linkage: WITH_LIBR=$WITH_LIBR =="

## BLAS/LAPACK, following the upstream recipe for each architecture:
## Apple Silicon uses ARMPL, Intel uses the Accelerate framework (which
## has its own code paths in gmrflib/simd.c). Either falls back to R's own
## Rblas/Rlapack when the preferred backend is unavailable.
STATIC_BLAS=${STATIC_BLAS:-1}
ARMPL_DIR=$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)
if [ "$ARCH" = arm64 ] && [ -n "$ARMPL_DIR" ]; then
    FLAGS="$FLAGS -DINLA_WITH_ARMPL -I$ARMPL_DIR/include"
    ## ARMPL's macOS build is flang-based: its static Fortran objects
    ## reference the flang runtime, which a GCC/gfortran link cannot supply
    ## (ilaenv_ and friends stay undefined), and they are compiled for macOS
    ## 13.0, which would silently raise this artifact's 11.0 floor. The
    ## dylib carries its own runtime, so ARMPL is linked dynamically here and
    ## ci/package-macos.sh copies it into the bundle: self-contained either
    ## way. On Linux the _gcc build is used and does embed (STATIC_BLAS=1).
    ARMPL_A=""
    [ "${ARMPL_MACOS_STATIC:-0}" = 1 ] && \
        ARMPL_A=$(ls "$ARMPL_DIR"/lib/libarmpl_lp64.a "$ARMPL_DIR"/lib/libarmpl.a 2>/dev/null | head -1 || true)
    if [ "$STATIC_BLAS" = 1 ] && [ -n "$ARMPL_A" ]; then
        echo "BLAS: ARMPL (static, embedded) $ARMPL_A"
        AMATH=$(ls "$ARMPL_DIR"/lib/libamath.a 2>/dev/null | head -1 || echo -lamath)
        ASTR=$(ls "$ARMPL_DIR"/lib/libastring.a 2>/dev/null | head -1 || echo -lastring)
        ## ARMPL's archive carries C++ built against libc++ (the std::__1
        ## symbols), while this binary is linked by g++ against libstdc++.
        ## Both can coexist -- they use different inline namespaces and
        ## ARMPL's C++ never crosses the boundary -- but libc++ has to be on
        ## the link line or those symbols stay undefined.
        BLASLIBS="$ARMPL_A $AMATH $ASTR -L$ARMPL_DIR/lib -lc++"
    else
        echo "BLAS: ARMPL (shared) at $ARMPL_DIR"
        BLASLIBS="-L$ARMPL_DIR/lib -Wl,-rpath,$ARMPL_DIR/lib -larmpl -lamath -lastring"
    fi
elif [ "$ARCH" = x86_64 ]; then
    echo "BLAS: Accelerate framework"
    FLAGS="$FLAGS -DINLA_WITH_FRAMEWORK_ACCELERATE"
    BLASLIBS="-framework Accelerate"
else
    echo "BLAS: R's Rblas/Rlapack"
    BLASLIBS="-lRblas -lRlapack"
fi

## ---- 1. External model packages --------------------------------------------
for d in "$EPATH"/*/; do
    [ -f "$d/Makefile" ] && make -C "$d" clean >/dev/null
done
rm -f "$EPATH"/lib*.a
( cd "$EPATH" && ./build \
      CC="$CC" CXX="$CXX" FC="$FC" \
      FLAGS="" \
      INC="-DINLA_WITH_EXTERNAL_PACKAGES -I$BREW/include/eigen3 \
           -I$RHOME/include -I$EPATH -I$ROOT/inlaprog/src" )

FAILED=0
for d in "$EPATH"/*/; do
    p=$(basename "$d")
    [ -f "$d/Makefile" ] || continue
    if [ ! -d "$d/$p" ]; then
        echo "WARNING: external package $p: repository not available, skipped"
    elif [ ! -f "$EPATH/lib$p.a" ]; then
        echo "ERROR: external package $p cloned but produced no archive"
        FAILED=1
    fi
done
[ "$FAILED" -eq 0 ]
ls "$EPATH"/lib*.a >/dev/null 2>&1 \
    || { echo "ERROR: no external-package archives at all (git access problem?)"; exit 1; }

## METIS 5.2 split its support routines into GKlib; link it when the build
## produced a separate archive.
GKLIB=""
[ -f "$DEPS/lib/libGKlib.a" ] && GKLIB="-lGKlib"

## ld64 whole-archive: one -force_load per archive.
EXTOBJ=""
for a in "$EPATH"/lib*.a; do
    EXTOBJ="$EXTOBJ -Wl,-force_load,$a"
done

## ---- 2. GMRFLib (+ vendored taucs/amd) --------------------------------------
## BSD ar has no 'U' flag; plain rv is deterministic enough here. SED=gsed:
## the dependency-file rule is written for GNU sed.
SED=$(command -v gsed || command -v sed)
make -C "$ROOT/gmrflib" -j"$JOBS" PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" ARFLAGS=rv SED="$SED"
make -C "$ROOT/gmrflib"           PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" ARFLAGS=rv SED="$SED" install

## ---- 3. inlaprog -> inla -----------------------------------------------------
## EXTLIBS3 is replaced below: its default is GNU-ld syntax
## (--whole-archive -lpthread), which ld64 rejects outright. Threads come
## from libSystem here, so nothing is lost.
make -C "$ROOT/inlaprog" -j"$JOBS" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" \
     EXTLIBS2="$EXTOBJ \
               -L$BREW/lib -L$DEPS/lib -L$BREW/opt/openssl@3/lib -L$BREW/opt/libtool/lib \
               -Wl,-rpath,$RHOME/lib -Wl,-rpath,$BREW/lib \
               -lgsl -lgslcblas $BLASLIBS \
               -lmuparser -lz -lmetis $GKLIB \
               -lltdl -lcrypto -lgfortran -lquadmath \
               -static-libstdc++ -static-libgcc -lm -ldl" \
     EXTLIBS3="-lm" \
     inla
make -C "$ROOT/inlaprog" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" install

## ---- 4. Runtime invariants ---------------------------------------------------
## One OpenMP runtime only. On macOS the danger is a mix of GCC's libgomp
## and LLVM's libomp: two runtimes in one process oversubscribe the machine
## and can abort. Keeping this single is also what allows the binary to be
## linked later against a library carrying its own BLAS.
OMPRT=$(otool -L "$PREFIX/bin/inla" 2>/dev/null | grep -oE 'lib(gomp|omp|iomp5)[^ ]*dylib' | sort -u || true)
echo "== OpenMP runtime: ${OMPRT:-none (static)} =="
if [ "$(echo "$OMPRT" | grep -c .)" -gt 1 ]; then
    echo "ERROR: more than one OpenMP runtime is linked:"; echo "$OMPRT"
    exit 1
fi

## ---- 5. Sanity --------------------------------------------------------------
"$PREFIX/bin/inla" -ping
echo "OK: inla built and installed in $PREFIX/bin (macOS arm64)"
