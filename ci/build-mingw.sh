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

## sTiles: an alternative sparse-matrix backend the sources already support
## (gmrflib/smtp-stiles.c, guarded by INLA_WITH_STILES). Only a prebuilt
## library and its header are needed -- ci/fetch-stiles.sh stages them from
## a published release, so no sTiles source is involved. It has to be
## resolved here, before GMRFLib is built: smtp-stiles.c is part of GMRFLib.
##
## The Windows asset ships an import library (libstiles.dll.a) beside the
## DLL, so this links the ordinary way. The DLL itself, and the runtime it
## brings, are picked up by the bundling step further down.
WITH_STILES=${WITH_STILES:-0}
STILES_LIBS=""
STILES_DIR=${STILES_DIR:-$PREFIX/stiles}
if [ "$WITH_STILES" = 1 ]; then
    [ -f "$STILES_DIR/include/stiles.h" ] \
        || { echo "ERROR: no stiles.h under $STILES_DIR (run ci/fetch-stiles.sh)"; exit 1; }
    IMP=$(ls "$STILES_DIR"/lib/libstiles.dll.a 2>/dev/null | head -1 || true)
    [ -n "$IMP" ] \
        || { echo "ERROR: the Windows libstiles asset has no import library (libstiles.dll.a)"; exit 1; }
    echo "== sTiles: $STILES_DIR =="
    ## -ltileindexer is NOT needed: released libraries carry it inside.
    FLAGS="$FLAGS -DINLA_WITH_STILES -I$STILES_DIR/include"
    STILES_LIBS="$IMP"
fi

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
               $LTDL $DL $STILES_LIBS \
               -lgfortran -lquadmath -lcrypto -lz \
               -lpthread -lm" \
     EXTLIBS3="-lm" \
     inla
[ -f "$ROOT/inlaprog/inla.exe" ] || { echo "ERROR: inla.exe was not produced"; exit 1; }

## No 32-bit pseudo-relocations may survive in the exe. They are created
## when code references DATA from a DLL without dllimport: the linker
## leaves a 32-bit slot the mingw runtime patches at startup, and with
## today's high-entropy ASLR the target lands beyond +-2GB often enough
## that the exe dies with "32 bit pseudo relocation ... out of range"
## before main() -- which is exactly how the first sTiles-linked build
## failed on a real Windows runner. --enable-auto-import was dropped from
## the link above so ld now WARNS per auto-imported symbol; this check
## makes the warnings fatal and names the culprits. (Read before strip:
## the list symbols are gone afterwards.)
PSTART=$($TRIPLET-nm "$ROOT/inlaprog/inla.exe" 2>/dev/null | awk '/ __RUNTIME_PSEUDO_RELOC_LIST__$/{print $1}')
PEND=$($TRIPLET-nm "$ROOT/inlaprog/inla.exe" 2>/dev/null | awk '/ __RUNTIME_PSEUDO_RELOC_LIST_END__$/{print $1}')
if [ -n "$PSTART" ] && [ "$PSTART" != "$PEND" ]; then
    echo "ERROR: the exe carries runtime pseudo-relocations (data auto-imports)."
    echo "       These crash at startup whenever ASLR places the DLL out of"
    echo "       32-bit range. The symbols involved (from the linker):"
    grep -iE "auto-import|pseudo" "$ROOT"/inlaprog/*.log 2>/dev/null || true
    exit 1
fi
echo "== pseudo-relocation check: clean =="

$TRIPLET-strip "$ROOT/inlaprog/inla.exe" 2>/dev/null || mingw-strip "$ROOT/inlaprog/inla.exe" || true
cp -f "$ROOT/inlaprog/inla.exe" "$PREFIX/bin/"

## ---- 4. Bundle: exe + every non-system DLL it imports ----------------------
OUT=$ROOT/dist-win
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$PREFIX/bin/inla.exe" "$OUT/"

## R's BLAS and LAPACK travel with the bundle. The exe imports them by name,
## so without them beside it Windows refuses to start the process at all --
## no message, just STATUS_DLL_NOT_FOUND (0xc0000135). Relying on the user's
## R to supply them only works when the exe is launched from R with R's bin
## directory on PATH; anyone unpacking the artifact and running inla.exe
## directly gets an exe that cannot start.
##
## Their dependencies come too, and that is not only these two files:
## Rblas.dll imports R.dll, which imports Rgraphapp.dll and Riconv.dll. The
## import walk below follows that chain, so the bundle ends up with the
## whole closure rather than the two names the exe happens to mention.
##
## Note what this does NOT change: the binary still has no static import of
## R.dll (that is INLA_WITH_LIBR_DLOPEN doing its job). R.dll is here to
## satisfy Rblas, and rgeneric still loads R at runtime from whichever R is
## driving the process.
if [ "$WITH_STILES" = 1 ]; then
    cp -v "$STILES_DIR"/lib/libstiles.dll "$OUT/"
fi

for dll in Rblas.dll Rlapack.dll; do
    [ -f "$RWIN/bin/x64/$dll" ] || { echo "ERROR: $dll not in $RWIN/bin/x64"; exit 1; }
    cp -v "$RWIN/bin/x64/$dll" "$OUT/"
done

## Is this name a Windows system DLL? Rather than guess from a list of
## names, ask the toolchain: mingw ships an import library (libfoo.a /
## libfoo.dll.a) for every DLL Windows itself provides, and for nothing
## else. R.dll alone imports a dozen of them (GDI32, COMDLG32, WINSPOOL,
## IMM32, ...), which is more than a hand-written list stays right about.
LIBDIRS="$SYSROOT/mingw/lib $SYSROOT2/mingw/lib /usr/x86_64-w64-mingw32/sys-root/mingw/lib /usr/x86_64-w64-mingw32/lib"
is_system_dll() {
    ## strip whatever extension it carries: R.dll imports WINSPOOL.DRV, and
    ## the import library for that one is still plain libwinspool.a
    b=$(echo "${1%.*}" | tr 'A-Z' 'a-z')
    case "$1" in api-ms-*|API-MS-*) return 0 ;; esac
    for d in $LIBDIRS; do
        [ -f "$d/lib$b.a" ] && return 0
        [ -f "$d/lib$b.dll.a" ] && return 0
    done
    return 1
}

## Resolve DLL imports from the sysroots, the deps tree and R's bin. Three
## passes: each one may add DLLs whose own imports the next pass resolves
## (inla.exe -> Rblas.dll -> R.dll -> Rgraphapp.dll/Riconv.dll is four
## levels, and the loop stops adding when a pass copies nothing new).
for pass in 1 2 3 4; do
    for exe in "$OUT"/*.exe "$OUT"/*.dll; do
        [ -f "$exe" ] || continue
        $TRIPLET-objdump -p "$exe" 2>/dev/null | awk '/DLL Name/ {print $3}'
    done | sort -u | while read -r dll; do
        is_system_dll "$dll" && continue
        [ -f "$OUT/$dll" ] && continue
        ## Both mingw sysroots: ltdl comes from the msvcrt one when the
        ## UCRT variant is not packaged (upstream links it the same way).
        ## $RWIN/bin/x64 supplies R's own DLLs -- Rblas and Rlapack, and
        ## through them R.dll, Rgraphapp.dll and Riconv.dll.
        ## $STILES_DIR/lib comes FIRST on purpose. libstiles.dll and this
        ## build share several runtime DLLs by name (libgomp-1.dll above
        ## all), and Windows loads one file per name from the exe's
        ## directory -- so whichever copy is staged is the copy BOTH use.
        ## sTiles is built with a newer GCC than the cross toolchain here,
        ## and libgomp is backward compatible, so its copy is the one that
        ## satisfies both. The reverse order can leave libstiles.dll unable
        ## to resolve a symbol its own build produced.
        src=$(find "$STILES_DIR/lib" "$SYSROOT/mingw/bin" "$SYSROOT2/mingw/bin" \
                   /usr/x86_64-w64-mingw32*/sys-root/mingw/bin \
                   /usr/lib/gcc/$TRIPLET "$DEPS" "$RWIN/bin/x64" \
                   -name "$dll" 2>/dev/null | head -1)
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
