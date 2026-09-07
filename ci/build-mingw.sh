#!/usr/bin/env bash
## Cross-compile inla.exe for Windows x86-64 with the MinGW UCRT toolchain,
## following the same stages as ci/build.sh: external packages -> GMRFLib ->
## inlaprog. Environment comes from ci/deps-mingw.sh. BLAS and LAPACK are a
## static OpenBLAS built by that script, deliberately NOT R's Rblas/Rlapack
## (the note there says why). In the default dlopen mode the exe references
## no R library at all, and rgeneric loads the running machine's R on first
## use, so exactly one R.dll is ever mapped into the process.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
## Versions are pinned in one file the maintainer owns; BUILDINFO below
## reports the BLAS from it rather than from a second copy of the number.
[ -f "$ROOT/ci/toolchain.env" ] && . "$ROOT/ci/toolchain.env"
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

## Version reported by the built binary. It is the R package's Version from
## rinla/DESCRIPTION, so `inla -V` and packageVersion("INLA") agree for anyone
## who installs the R package and the binary from the same commit. The short
## commit is kept in INLA_TAG as <version>+<sha>. No spaces or parens:
## this string travels through a $(MAKE) recipe into /bin/sh, where dash
## treats a bare "(" as a syntax error and the whole build dies.
## which is where the build-traceability belongs; GITCOMMIT has to stay a bare
## preprocessor token, so it carries the version alone.
SHA=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
TAG=$(sed -n 's/^Version:[[:space:]]*//p' "$ROOT/rinla/DESCRIPTION" 2>/dev/null | head -1)
## Reported VERBATIM from DESCRIPTION, zeros and all: "26.09.03", matching the
## release tag (v26.09.03), the startup banner and DESCRIPTION itself. Only
## packageVersion("INLA") differs, showing "26.9.3", because R's numeric_version
## strips leading zeros and cannot be told otherwise. Comparisons are unaffected:
## the version check uses package_version(), which normalises both sides.
[ -n "$TAG" ] || TAG=$SHA
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
 -DGITCOMMIT=$TAG -DINLA_TAG='\"$TAG+$SHA\"' \
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
## The DLL and the runtime it brings are picked up by the bundling step
## further down.
WITH_STILES=${WITH_STILES:-0}
STILES_LIBS=""
STILES_DIR=${STILES_DIR:-$PREFIX/stiles}
if [ "$WITH_STILES" = 1 ]; then
    [ -f "$STILES_DIR/include/stiles.h" ] \
        || { echo "ERROR: no stiles.h under $STILES_DIR (run ci/fetch-stiles.sh)"; exit 1; }
    ## Link against the DLL ITSELF rather than its import library. The
    ## shipped import library (.dll.a) types every
    ## function as DATA (a PE --version-script on the DLL link breaks
    ## --out-implib's classification), so using it leaves the exe full of
    ## 32-bit runtime pseudo-relocations that crash at startup under ASLR.
    ## ld reading the DLL's own export table sees the symbols in .text and
    ## synthesizes correct thunks -- verified: the pseudo-reloc list is
    ## empty, with and without LTO.
    [ -f "$STILES_DIR/lib/libstiles.dll" ] \
        || { echo "ERROR: no libstiles.dll under $STILES_DIR/lib"; exit 1; }
    ## -ltileindexer is NOT needed: released libraries carry it inside.
    FLAGS="$FLAGS -DINLA_WITH_STILES -I$STILES_DIR/include"
    STILES_LIBS="$STILES_DIR/lib/libstiles.dll"
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
## static OpenBLAS/Rmath/gsl/metis/muparser,
## static C++/GCC runtimes, ltdl import library from the sysroot.
## The GCC runtime is NOT linked with a blanket -static: libgomp must stay
## a DLL so that a library loaded into this process later (libstiles.dll,
## which is built the same way) SHARES one OpenMP runtime with the binary.
## Two copies of libgomp in one process is the "OMP: Error #15" case. The
## bundling step below ships libgomp-1.dll and friends beside the exe, and
## now fails loudly if it cannot find one.
## No libgslcblas: OpenBLAS exports the cblas_* interface GSL calls, exactly
## as MKL and Accelerate do on the other platforms. It was needed only while
## the BLAS was R's Rblas, which is Fortran-only.
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

## mimalloc first on the link line + --undefined=mi_version, exactly the
## upstream Windows recipe: forcing mi_version in pulls mimalloc's malloc
## override ahead of the UCRT one.
MIMALLOC=""
[ -f "$DEPS/lib/libmimalloc.dll.a" ] && MIMALLOC="$DEPS/lib/libmimalloc.dll.a -Wl,--undefined=mi_version"

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
     EXTLIBS2="$MIMALLOC \
               -Wl,--whole-archive $EXTOBJ -Wl,--no-whole-archive \
               -static-libstdc++ -static-libgcc \
               $DEPS/lib/libRmath.a $DEPS/lib/libgsl.a \
               $DEPS/lib/libmetis.a $GKLIB $DEPS/lib/libmuparser.dll.a \
               $DEPS/lib/libopenblas.a \
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
    ## The workflow runs this script with all output redirected to
    ## $ROOT/build.log, so the linker's per-symbol auto-import warnings are
    ## in that file by the time the link has finished.
    grep -iE "auto-import|pseudo-reloc" "$ROOT/build.log" 2>/dev/null | head -30 || true
    exit 1
fi
echo "== pseudo-relocation check: clean =="

$TRIPLET-strip "$ROOT/inlaprog/inla.exe" 2>/dev/null || mingw-strip "$ROOT/inlaprog/inla.exe" || true
cp -f "$ROOT/inlaprog/inla.exe" "$PREFIX/bin/"

## ---- 4. Bundle: exe + every non-system DLL it imports ----------------------
OUT=$ROOT/dist-win
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$PREFIX/bin/inla.exe" "$OUT/"

## No R library travels with the bundle in dlopen mode, and that is the
## whole point. The exe used to import R's Rblas/Rlapack, which import
## R.dll, so the bundle had to ship R.dll beside the exe; Windows resolves a
## dependency from the executable's directory before PATH, so that copy was
## mapped at process start, and rgeneric then dlopened the user's R.dll from
## R_HOME. Two R.dll copies in one process segfault the embedded R as soon
## as an rgeneric model is fitted, even at identical versions. With a static
## OpenBLAS the exe references no R library, so the only R.dll in the
## process is the one rgeneric loads. The guard after the import walk keeps
## it that way.
if [ "$WITH_STILES" = 1 ]; then
    cp -v "$STILES_DIR"/lib/libstiles.dll "$OUT/"
fi

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

## Resolve DLL imports from the sysroots and the deps tree. Several passes:
## each one may add DLLs whose own imports the next pass resolves, and the
## loop stops adding when a pass copies nothing new.
##
## R's bin is searched ONLY when R.dll is linked at build time (WITH_LIBR=1).
## In dlopen mode it must not be, or the walk would quietly pull R's DLLs
## back into the bundle and reintroduce the two-R.dll crash.
RSRC=""
[ "$WITH_LIBR" = 2 ] || RSRC="$RWIN/bin/x64"
for pass in 1 2 3 4; do
    for exe in "$OUT"/*.exe "$OUT"/*.dll; do
        [ -f "$exe" ] || continue
        $TRIPLET-objdump -p "$exe" 2>/dev/null | awk '/DLL Name/ {print $3}'
    done | sort -u | while read -r dll; do
        is_system_dll "$dll" && continue
        [ -f "$OUT/$dll" ] && continue
        ## Both mingw sysroots: ltdl comes from the msvcrt one when the
        ## UCRT variant is not packaged (upstream links it the same way).
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
                   /usr/lib/gcc/$TRIPLET "$DEPS" $RSRC \
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

## The regression guard for the two-R.dll crash. In dlopen mode neither the
## exe nor anything beside it may reference an R library: if one does,
## Windows maps it at startup and rgeneric's own dlopen adds a second copy,
## which segfaults the embedded R. Checked here rather than trusted, because
## the failure is silent until someone fits an rgeneric model on Windows.
if [ "$WITH_LIBR" = 2 ]; then
    BAD=$($TRIPLET-objdump -p "$OUT"/*.exe "$OUT"/*.dll 2>/dev/null \
          | awk '/DLL Name/ {print $3}' \
          | grep -iE '^(R|Rblas|Rlapack|Rgraphapp|Riconv)\.dll$' | sort -u || true)
    if [ -n "$BAD" ]; then
        echo "ERROR: the bundle references an R library in dlopen mode:"
        printf '  %s\n' $BAD
        exit 1
    fi
    if ls "$OUT"/R*.dll >/dev/null 2>&1; then
        echo "ERROR: the bundle ships an R library in dlopen mode:"
        ls -1 "$OUT"/R*.dll
        exit 1
    fi
    echo "OK: no R library in the bundle; rgeneric will load exactly one R.dll"
fi

bash "$ROOT/ci/write-buildinfo.sh" "$OUT/BUILDINFO" "$CC" "$FLAGS" \
     "OpenBLAS ${OPENBLAS_VERSION:-0.3.29}, static, runtime CPU dispatch"
echo "== imports of the shipped exe =="
$TRIPLET-objdump -p "$OUT/inla.exe" | awk '/DLL Name/ {print "  " $3}' | sort -u
echo "== bundled files =="
ls -1 "$OUT"

( cd "$ROOT" && zip -qr inla-windows-x86_64.zip "$(basename "$OUT")" )
echo "OK: $(du -h "$ROOT/inla-windows-x86_64.zip" | cut -f1) Windows bundle"
