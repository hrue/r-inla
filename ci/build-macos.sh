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

FLAGS="$OPTFLAGS -mcpu=apple-m1 -pipe -pthread \
 -fopenmp -fopenmp-simd -flax-vector-conversions \
 -DINLA_WITH_SIMDE -DINLA_WITH_DEVEL \
 -DINLA_WITH_EXTERNAL_PACKAGES -DINLA_WITH_MUPARSER \
 -DGITCOMMIT=$TAG \
 -I$BREW/include -I$DEPS/include"

mkdir -p "$PREFIX"/bin "$PREFIX"/lib "$PREFIX"/include "$PREFIX/include.boot"
ln -sfn "$ROOT/gmrflib" "$PREFIX/include.boot/GMRFLib"
FLAGS="$FLAGS -I$PREFIX/include.boot"

RLIB_INC="-DINLA_WITH_LIBR -I$RHOME/include"
## -L$DEPS/lib supplies the standalone Rmath built by ci/deps-macos.sh: the
## framework ships libR but not that one.
RLIB_LIB="-L$RHOME/lib -L$DEPS/lib -lR -lRmath"

## BLAS/LAPACK: ARMPL when installed (what the upstream macOS binaries
## use; the binary then expects ARMPL on the running machine, like
## upstream's), otherwise R's own Rblas/Rlapack.
ARMPL_DIR=$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)
if [ -n "$ARMPL_DIR" ]; then
    echo "BLAS: ARMPL at $ARMPL_DIR"
    FLAGS="$FLAGS -DINLA_WITH_ARMPL -I$ARMPL_DIR/include"
    BLASLIBS="-L$ARMPL_DIR/lib -Wl,-rpath,$ARMPL_DIR/lib -larmpl -lamath -lastring"
else
    echo "BLAS: R's Rblas/Rlapack (ARMPL not found)"
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
               -lmuparser -lz -lmetis \
               -lltdl -lcrypto -lgfortran -lquadmath \
               -static-libstdc++ -static-libgcc -lm -ldl" \
     EXTLIBS3="-lm" \
     inla
make -C "$ROOT/inlaprog" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" install

## ---- 4. Sanity --------------------------------------------------------------
"$PREFIX/bin/inla" -ping
echo "OK: inla built and installed in $PREFIX/bin (macOS arm64)"
