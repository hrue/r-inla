#!/usr/bin/env bash
## Cross-compile inla.exe for Windows x86-64 with the MinGW UCRT toolchain,
## following the same stages as ci/build.sh: external packages -> GMRFLib ->
## inlaprog. Environment comes from ci/deps-mingw.sh. The BLAS/LAPACK are
## the ones Windows R itself ships (Rblas.dll/Rlapack.dll), and libR is the
## real R.dll, so the exe runs against the user's installed R for Windows.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local-win}
DEPS=${DEPS:-/opt/mingw-deps}
JOBS=${JOBS:-$(nproc)}

TRIPLET=x86_64-w64-mingw32ucrt
CC=$TRIPLET-gcc
CXX=$TRIPLET-g++
FC=$TRIPLET-gfortran
SYSROOT=$($CC --print-sysroot)
SYSROOT2=$(x86_64-w64-mingw32-gcc --print-sysroot 2>/dev/null || echo "$SYSROOT")
EPATH=$ROOT/external-packages
RWIN=$DEPS/R-win

TAG=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo devel)
echo "== building $TAG for Windows with $CC ($($CC --version | head -1)) =="

## No NUMA/CLONE_TARGETS (Linux-only); otherwise the devel feature set.
## Optimization: the upstream Windows configuration. LTO stays opt-in.
OPT=${OPT:-fast}
LTO=${LTO:-0}
case "$OPT" in
    fast) OPTFLAGS="-O3 -ftree-vectorize -funroll-loops -fvariable-expansion-in-unroller -ftracer" ;;
    safe) OPTFLAGS="-O2 -ftree-vectorize" ;;
    *)    echo "ERROR: OPT must be fast or safe"; exit 1 ;;
esac
[ "$LTO" = 1 ] && OPTFLAGS="$OPTFLAGS -flto=auto -ffat-lto-objects"
echo "== optimization: OPT=$OPT LTO=$LTO =="

FLAGS="$OPTFLAGS -mtune=generic -pipe -pthread -Wall -Wextra \
 -fopenmp -fopenmp-simd -flax-vector-conversions \
 -DINLA_WITH_SIMDE -DINLA_WITH_DEVEL \
 -DINLA_WITH_EXTERNAL_PACKAGES -DINLA_WITH_MUPARSER \
 -DGITCOMMIT=$TAG -DINLA_TAG='\"$TAG\"' \
 -I$DEPS/include -I$RWIN/include"

mkdir -p "$PREFIX"/bin "$PREFIX"/lib "$PREFIX"/include
mkdir -p "$PREFIX/include.boot"
ln -sfn "$ROOT/gmrflib" "$PREFIX/include.boot/GMRFLib"
FLAGS="$FLAGS -I$PREFIX/include.boot"

## ---- 1. External model packages --------------------------------------------
for d in "$EPATH"/*/; do
    [ -f "$d/Makefile" ] && make -C "$d" clean >/dev/null
done
rm -f "$EPATH"/lib*.a
( cd "$EPATH" && ./build \
      CC="$CC" CXX="$CXX" FC="$FC" \
      FLAGS="" \
      INC="-DINLA_WITH_EXTERNAL_PACKAGES -I$DEPS/include -I$DEPS/include/eigen3 \
           -I$SYSROOT2/mingw/include/eigen3 -I$SYSROOT/mingw/include/eigen3 \
           -I$RWIN/include -I$EPATH -I$ROOT/inlaprog/src" )

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

## ---- 2. GMRFLib (+ vendored taucs/amd) --------------------------------------
make -C "$ROOT/gmrflib" -j"$JOBS" PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" AR=$TRIPLET-ar ARFLAGS=rvU
make -C "$ROOT/gmrflib"           PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" AR=$TRIPLET-ar install

## ---- 3. inlaprog -> inla.exe -------------------------------------------------
## Link model follows the upstream Windows recipe: whole-archive externals,
## Windows R's own R.dll/Rblas/Rlapack, static Rmath/gsl/metis/muparser,
## static C++/GCC runtimes, ltdl import library from the sysroot.
## The GCC runtime is NOT linked with a blanket -static: libgomp must stay
## a DLL so that a library loaded into this process later (libstiles.dll,
## which is built the same way) SHARES one OpenMP runtime with the binary.
## Two copies of libgomp in one process is the "OMP: Error #15" case. The
## bundling step below ships libgomp-1.dll and friends beside the exe, and
## now fails loudly if it cannot find one.
## libgslcblas is needed here although Linux and macOS get by without it:
## OpenBLAS, MKL and Accelerate all export the cblas_* interface GSL calls,
## while R's Rblas is Fortran-only.
## R linkage: 1 links R.dll at build time (what upstream ships), 2 loads
## the running machine's R through libltdl on first rgeneric use.
WITH_LIBR=${WITH_LIBR:-1}
if [ "$WITH_LIBR" = 2 ]; then
    ## The R include path stays: R.dll is not linked, but rmath.h still
    ## includes <Rmath.h> for the standalone math library.
    RLIB_INC="-DINLA_WITH_LIBR -DINLA_WITH_LIBR_DLOPEN -I$RWIN/include"
    RLIB_LIB=""
else
    RLIB_INC="-DINLA_WITH_LIBR -I$RWIN/include"
    RLIB_LIB="$RWIN/bin/x64/R.dll"
fi
echo "== R linkage: WITH_LIBR=$WITH_LIBR =="

## METIS 5.2 split its support routines into GKlib; it must follow libmetis
## on the link line.
GKLIB=""
[ -f "$DEPS/lib/libGKlib.a" ] && GKLIB="$DEPS/lib/libGKlib.a"

EXTOBJ=$(echo "$EPATH"/lib*.a)
LTDL=$(ls "$DEPS"/lib/libltdl.dll.a "$SYSROOT"/mingw/lib/libltdl.dll.a \
          "$SYSROOT2"/mingw/lib/libltdl.dll.a 2>/dev/null | head -1 || true)
## dlfcn-win32: from the sysroot when packaged, else the static build in $DEPS
DL=$(ls "$SYSROOT"/mingw/lib/libdl.dll.a "$SYSROOT2"/mingw/lib/libdl.dll.a \
        "$DEPS"/lib/libdl.a 2>/dev/null | head -1 || true)
make -C "$ROOT/inlaprog" -j"$JOBS" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" \
     RLIB_LIB="$RLIB_LIB" \
     EXTLIBS2="-Wl,--whole-archive $EXTOBJ -Wl,--no-whole-archive \
               -static-libstdc++ -static-libgcc \
               $DEPS/lib/libRmath.a $DEPS/lib/libgsl.a $DEPS/lib/libgslcblas.a \
               $DEPS/lib/libmetis.a $GKLIB $DEPS/lib/libmuparser.dll.a \
               $RWIN/bin/x64/Rblas.dll $RWIN/bin/x64/Rlapack.dll \
               $LTDL $DL \
               -lgfortran -lquadmath -lcrypto -lz \
               -Wl,--enable-auto-import -lpthread -lm" \
     EXTLIBS3="-lm" \
     inla
[ -f "$ROOT/inlaprog/inla.exe" ] || { echo "ERROR: inla.exe was not produced"; exit 1; }
$TRIPLET-strip "$ROOT/inlaprog/inla.exe" 2>/dev/null || mingw-strip "$ROOT/inlaprog/inla.exe" || true
cp -f "$ROOT/inlaprog/inla.exe" "$PREFIX/bin/"

## ---- 4. Bundle: exe + every non-system DLL it imports ----------------------
OUT=$ROOT/dist-win
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$PREFIX/bin/inla.exe" "$OUT/"
## Resolve DLL imports from the sysroots and the deps tree; R's own DLLs
## (R.dll, Rblas, Rlapack) are deliberately NOT bundled: they come from the
## user's installed R for Windows, found via R's bin directory on PATH.
for pass in 1 2 3; do
    for exe in "$OUT"/*.exe "$OUT"/*.dll; do
        [ -f "$exe" ] || continue
        $TRIPLET-objdump -p "$exe" 2>/dev/null | awk '/DLL Name/ {print $3}'
    done | sort -u | while read -r dll; do
        case "$dll" in
            R.dll|Rblas.dll|Rlapack.dll|KERNEL32*|msvcrt*|api-ms-*|ucrtbase*|ADVAPI32*|WS2_32*|USER32*|SHELL32*|ole32*|RPCRT4*|dbghelp*|bcrypt*|CRYPT32*) continue ;;
        esac
        [ -f "$OUT/$dll" ] && continue
        src=$(find "$SYSROOT/mingw/bin" "$SYSROOT2/mingw/bin" \
                   /usr/lib/gcc/$TRIPLET "$DEPS" -name "$dll" 2>/dev/null | head -1)
        if [ -n "$src" ]; then
            cp -v "$src" "$OUT/"
        else
            ## Never skip quietly: a missing DLL makes the exe fail to start
            ## on Windows with a bare exit code and no message at all.
            echo "MISSING: $dll (imported but not found on this system)" >> "$OUT/.missing"
        fi
    done
done

if [ -s "$OUT/.missing" ]; then
    echo "ERROR: the bundle is incomplete:"
    sort -u "$OUT/.missing"
    rm -f "$OUT/.missing"
    exit 1
fi
rm -f "$OUT/.missing"

echo "== imports of the shipped exe =="
$TRIPLET-objdump -p "$OUT/inla.exe" | awk '/DLL Name/ {print "  " $3}' | sort -u
echo "== bundled files =="
ls -1 "$OUT"

( cd "$ROOT" && zip -qr inla-windows-x86_64.zip "$(basename "$OUT")" )
echo "OK: $(du -h "$ROOT/inla-windows-x86_64.zip" | cut -f1) Windows bundle"
