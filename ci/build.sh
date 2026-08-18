#!/usr/bin/env bash
## Build the inla binary from this source tree, using only the makefiles
## already present (external-packages/build, gmrflib/Makefile,
## inlaprog/Makefile). All configuration is passed as make command-line
## overrides; no makefile is edited.
##
## Usage:  bash ci/build.sh          (system deps: see ci/deps-ubuntu.sh)
##
## Knobs (environment variables):
##   PREFIX      install prefix                     default: <repo>/local
##   WITH_LIBR   1 = embed R, enables rgeneric      default: 1
##   JOBS        parallel compile jobs              default: nproc
##   CC/CXX/FC   compilers                          default: gcc/g++/gfortran
##   BLAS        openblas | armpl                   default: openblas
##   OPT         fast | safe                        default: fast
##   LTO         1 to add -flto=auto                default: 0
##   WITH_STILES 1 to link libstiles (see ci/fetch-stiles.sh)  default: 0
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
WITH_LIBR=${WITH_LIBR:-1}
JOBS=${JOBS:-$(nproc)}
## Compilers: the newest versioned GCC installed, else the default one.
## (A container that puts its own gcc first -- the manylinux gcc-toolset,
## for instance -- keeps winning, since no gcc-N binaries exist there.)
newest_gcc() {
    ls /usr/bin/"$1"-1[0-9] 2>/dev/null | sort -V | tail -1 || true
}
CC=${CC:-$(newest_gcc gcc)}
CXX=${CXX:-$(newest_gcc g++)}
FC=${FC:-$(newest_gcc gfortran)}
CC=${CC:-gcc}
CXX=${CXX:-g++}
FC=${FC:-gfortran}
BLAS=${BLAS:-openblas}
EPATH=$ROOT/external-packages

## Optimization. "fast" is the upstream release configuration; "safe"
## exists for bisecting a problem that only appears optimized. LTO is the
## one upstream flag kept opt-in: it roughly doubles link time and its
## failures are unrelated to the change being tested.
OPT=${OPT:-fast}
LTO=${LTO:-0}
case "$OPT" in
    fast) OPTFLAGS="-O3 -ftree-vectorize -funroll-loops -fvariable-expansion-in-unroller -ftracer" ;;
    safe) OPTFLAGS="-O2 -ftree-vectorize" ;;
    *)    echo "ERROR: OPT must be fast or safe"; exit 1 ;;
esac
[ "$LTO" = 1 ] && OPTFLAGS="$OPTFLAGS -flto=auto -ffat-lto-objects"

## Version string compiled into the binary (shown by inla -V and inla -ping);
## falls back to "devel" outside a git checkout.
TAG=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo devel)

echo "== building $TAG with CC=$CC CXX=$CXX FC=$FC ($($CC --version | head -1)) =="
echo "== optimization: OPT=$OPT LTO=$LTO =="

## BLAS/LAPACK. Linked STATICALLY wherever an archive exists, so the BLAS
## ends up inside the binary as it does upstream, instead of travelling
## beside it. Set STATIC_BLAS=0 to link dynamically instead.
##
## The defines are not cosmetic: the sources have dedicated code paths
## behind INLA_WITH_MKL and INLA_WITH_ARMPL (dot products, idxval,
## lapack-interface), so each belongs with its libraries.
STATIC_BLAS=${STATIC_BLAS:-1}
case "$BLAS" in
    mkl)
        MKLROOT=${MKLROOT:-/opt/intel/oneapi/mkl/latest}
        MKLLIB=$MKLROOT/lib/intel64
        [ -d "$MKLLIB" ] || MKLLIB=$MKLROOT/lib
        [ -f "$MKLLIB/libmkl_core.a" ] || { echo "ERROR: no static MKL under $MKLLIB"; exit 1; }
        echo "== BLAS: MKL (static, embedded) from $MKLLIB =="
        BLAS_INC="-DINLA_WITH_MKL -I$MKLROOT/include"
        ## The start/end group resolves MKL's circular dependencies, as in
        ## the upstream recipe.
        BLAS_LIBS="-Wl,--start-group $MKLLIB/libmkl_intel_lp64.a \
                   $MKLLIB/libmkl_gnu_thread.a $MKLLIB/libmkl_core.a -Wl,--end-group"
        ;;
    armpl)
        ARMPL_DIR=${ARMPL_DIR:-$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)}
        [ -n "$ARMPL_DIR" ] || { echo "ERROR: BLAS=armpl but no /opt/arm/armpl_* found"; exit 1; }
        BLAS_INC="-DINLA_WITH_ARMPL -I$ARMPL_DIR/include"
        ARMPL_A=$(ls "$ARMPL_DIR"/lib/libarmpl_lp64.a "$ARMPL_DIR"/lib/libarmpl.a 2>/dev/null | head -1 || true)
        if [ "$STATIC_BLAS" = 1 ] && [ -n "$ARMPL_A" ]; then
            echo "== BLAS: ARMPL (static, embedded) $ARMPL_A =="
            AMATH=$(ls "$ARMPL_DIR"/lib/libamath.a 2>/dev/null | head -1 || echo -lamath)
            ## libastring is deliberately NOT linked. Its vectorized memchr
            ## reads 16 bytes at a time, and ARMPL's own /proc/cpuinfo parser
            ## (a static initializer, before main) calls it on a small heap
            ## buffer: when that buffer lies near a page end the over-read
            ## touches an unmapped page and the process dies. Measured at
            ## ~3% of startups on the arm64 runners, in every ARMPL version
            ## tested (26.07, 25.07, 25.04); linking libarmpl WITHOUT
            ## libastring is clean under valgrind. Only string routines are
            ## lost, which INLA does not use in any hot path.
            BLAS_LIBS="$ARMPL_A $AMATH -L$ARMPL_DIR/lib"
        else
            echo "== BLAS: ARMPL (shared) at $ARMPL_DIR =="
            BLAS_LIBS="-L$ARMPL_DIR/lib -larmpl -lamath"
        fi
        ;;
    *)
        BLAS_INC=""
        OB_A=$(ls /usr/lib64/libopenblas.a /usr/lib/*-linux-gnu/libopenblas.a 2>/dev/null | head -1 || true)
        if [ "$STATIC_BLAS" = 1 ] && [ -n "$OB_A" ]; then
            echo "== BLAS: OpenBLAS (static, embedded) $OB_A =="
            BLAS_LIBS="$OB_A"
        else
            echo "== BLAS: OpenBLAS (shared) =="
            BLAS_LIBS="-lopenblas"
        fi
        ;;
esac

## sTiles: an alternative sparse-matrix backend the sources already support
## (gmrflib/smtp-stiles.c, guarded by INLA_WITH_STILES). Only a prebuilt
## library and its header are needed -- ci/fetch-stiles.sh stages them from
## a published release, so no sTiles source is involved.
WITH_STILES=${WITH_STILES:-0}
STILES_LIBS=""
if [ "$WITH_STILES" = 1 ]; then
    STILES_DIR=${STILES_DIR:-$PREFIX/stiles}
    [ -f "$STILES_DIR/include/stiles.h" ] \
        || { echo "ERROR: no stiles.h under $STILES_DIR (run ci/fetch-stiles.sh)"; exit 1; }
    echo "== sTiles: $STILES_DIR =="
    ## -ltileindexer is NOT needed: released libraries carry it inside.
    STILES_INC="-DINLA_WITH_STILES -I$STILES_DIR/include"
    STILES_LIBS="-L$STILES_DIR/lib -Wl,-rpath,$STILES_DIR/lib -lstiles"
else
    STILES_INC=""
fi

## The devel feature set, minus PARDISO (requires a license; libpardiso.c
## provides stubs so the link still closes) and minus MKL (this build links
## OpenBLAS). INLA_WITH_SIMDE needs the SIMDE headers (libsimde-dev).
## No -march: baseline x86-64, with the fast per-CPU paths supplied at run
## time by the INLA_WITH_CLONE_TARGETS function clones.
## MARCH: default is the universal x86-64 baseline (runtime dispatch via
## MKL + INLA_WITH_CLONE_TARGETS covers the hot paths). The modern lane
## passes -march=x86-64-v3 so ALL code gets AVX2/FMA, closing most of the
## portable-vs-native gap on any CPU from ~2015 on.
MARCH=${MARCH:--mtune=generic}
FLAGS="$OPTFLAGS $MARCH -pipe -pthread -Wall -Wextra \
 -fopenmp -fopenmp-simd -flax-vector-conversions \
 -DINLA_WITH_SIMDE -DINLA_WITH_DEVEL -DINLA_WITH_CLONE_TARGETS \
 -DINLA_WITH_EXTERNAL_PACKAGES -DINLA_WITH_MUPARSER -DINLA_WITH_NUMA \
 -DGITCOMMIT=$TAG -DINLA_TAG='\"$TAG\"' $BLAS_INC $STILES_INC"

## R linkage, three modes:
##   WITH_LIBR=1  link the shared libR at build time (rgeneric works;
##                needs R with a shared libR on the BUILD machine)
##   WITH_LIBR=2  dlopen mode: no libR at build or start time; rgeneric
##                loads the running machine's libR (via R_HOME) on first
##                use. Needs no R headers, so it builds without R.
##   WITH_LIBR=0  no R embedding at all; everything but rgeneric works.
if [ "$WITH_LIBR" = 1 ]; then
    ## libR's directory differs per family (Debian: /usr/lib/R/lib,
    ## Fedora: /usr/lib64/R/lib); ask R itself.
    R_HOME_DIR=${R_HOME_DIR:-$(R RHOME 2>/dev/null || echo /usr/lib/R)}
    RLIB_INC="-DINLA_WITH_LIBR -I/usr/include/R -I/usr/share/R/include"
    RLIB_LIB="-lRmath -L$R_HOME_DIR/lib -lR"
elif [ "$WITH_LIBR" = 2 ]; then
    ## The R include path stays: libR is not linked, but rmath.h still
    ## includes <Rmath.h> for the standalone math library.
    RLIB_INC="-DINLA_WITH_LIBR -DINLA_WITH_LIBR_DLOPEN -I/usr/include/R -I/usr/share/R/include"
    RLIB_LIB="-lRmath"
else
    RLIB_INC=""
    RLIB_LIB="-lRmath"
fi

mkdir -p "$PREFIX"/bin "$PREFIX"/lib "$PREFIX"/include

## Header bootstrap. gmrflib's headers reference each other as
## "GMRFLib/...", which on an installed system resolves through
## $PREFIX/include/GMRFLib -- i.e. through a PREVIOUS install. A fresh
## checkout has no install yet (gmrflib/fsort/fluxsort.h -> quadsort.h is
## the include that trips), so provide the same mapping via a private
## include dir holding one symlink to the source tree.
mkdir -p "$PREFIX/include.boot"
ln -sfn "$ROOT/gmrflib" "$PREFIX/include.boot/GMRFLib"
FLAGS="$FLAGS -I$PREFIX/include.boot"

## ---- 1. External model packages -------------------------------------------
## Clones each package (network), compiles it, and produces one lib<pkg>.a
## per package in external-packages/ plus the generated cgeneric-table.h/
## cgeneric-defs.h. Built in-tree; absolute include paths (the compile runs
## three directories deep inside each cloned package, so relative ones are
## fragile). The R include dirs cover both layouts (Fedora: /usr/include/R,
## Debian/Ubuntu: /usr/share/R/include) since the package sources use R.h
## and Rmath.h. -fPIC so the objects can also serve in a shared library.
##
## Start from clean package dirs: with a leftover clone, "make download"
## fails and the build silently reuses stale sources. Old archives go too,
## so a package that fails to build cannot be masked by a stale lib<pkg>.a
## from an earlier run.
for d in "$EPATH"/*/; do
    [ -f "$d/Makefile" ] && make -C "$d" clean >/dev/null
done
rm -f "$EPATH"/lib*.a
( cd "$EPATH" && ./build \
      CC="$CC" CXX="$CXX" FC="$FC" \
      FLAGS="-fPIC" \
      INC="-DINLA_WITH_EXTERNAL_PACKAGES \
           -I/usr/local/include/eigen3 -I/usr/include/eigen3 \
           -I/usr/include/R -I/usr/share/R/include \
           -I$EPATH -I$ROOT/inlaprog/src" )

## The build script above does not stop on a failed package, so check the
## outcome per package. A package whose repository cannot be cloned at all
## is skipped with a warning: its entries never reach the generated model
## table, so the link stays consistent without it. A package that cloned
## but produced no archive is a compile failure and stops the build.
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

## ---- 2. GMRFLib (+ vendored taucs/amd) -> $PREFIX/lib/libGMRFLib.a ---------
## ARFLAGS=rvU: modern ar defaults to deterministic mode; U restores real
## timestamps so make's staleness checks work across build and install.
## Note: this tree has no gmrflib/doc directory, so the install step prints
## a harmless rsync error for it and still succeeds (the recipe ends in echo).
make -C "$ROOT/gmrflib" -j"$JOBS" PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" ARFLAGS=rvU
make -C "$ROOT/gmrflib"           PREFIX="$PREFIX" FLAGS="$FLAGS" \
     CC="$CC" CXX="$CXX" FC="$FC" install

## ---- 3. inlaprog -> $PREFIX/bin/inla ---------------------------------------
## Link-line notes, relative to the Makefile's defaults:
##   - the six external archives go in with --whole-archive so their model
##     registration objects are not dropped by the linker
##   - openblas replaces "-llapack -lblas -lgslcblas": one BLAS+LAPACK+CBLAS
##     implementation in the process instead of three partial ones
##   - numa/hwloc (my-numa.c) and ltdl (cgeneric loading in inla.c) are used
##     by the sources but absent from the default link line
##   - -ldl stays although our own loading goes through ltdl: the external
##     packages call dl* directly, and glibc before 2.34 (e.g. the
##     manylinux container) ships libdl as a separate library
## -I$EPATH lets inla.c find the generated cgeneric-table.h.
## One space-separated line (ls would emit newlines, and a newline inside a
## make command-line variable truncates the link recipe at that point,
## leaving --whole-archive unterminated).
EXTOBJ=$(echo "$EPATH"/lib*.a)
## libquadmath is x86-only (__float128); aarch64 neither has nor needs it.
QUADMATH=""
[ "$(uname -m)" = "x86_64" ] && QUADMATH="-lquadmath"
make -C "$ROOT/inlaprog" -j"$JOBS" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" \
     EXTLIBS2="-Wl,--whole-archive $EXTOBJ -Wl,--no-whole-archive \
               $STILES_LIBS \
               -lgsl $BLAS_LIBS -lmuparser -lz -lmetis \
               -lnuma -lhwloc -lltdl -lcrypto -lgfortran $QUADMATH -lm -ldl"
make -C "$ROOT/inlaprog" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" install

## ---- 4. Runtime invariants ---------------------------------------------------
## Exactly one OpenMP runtime, and it must be libgomp. A process holding
## both libgomp and Intel's libiomp5 oversubscribes its threads and can
## crash outright, and that is precisely what happens if MKL is linked
## through its Intel-threaded variant instead of mkl_gnu_thread. It is also
## the invariant to keep when this binary is later linked against a library
## that carries its own BLAS.
OMPRT=$(ldd "$PREFIX/bin/inla" 2>/dev/null | grep -oE 'lib(gomp|iomp5|omp)\.so[^ ]*' | sort -u || true)
echo "== OpenMP runtime: ${OMPRT:-none (static)} =="
if echo "$OMPRT" | grep -q 'libiomp5'; then
    echo "ERROR: Intel's libiomp5 is linked; expected libgomp"
    ldd "$PREFIX/bin/inla" | grep -iE 'omp|mkl'
    exit 1
fi
if [ "$(echo "$OMPRT" | grep -c .)" -gt 1 ]; then
    echo "ERROR: more than one OpenMP runtime is linked:"; echo "$OMPRT"
    exit 1
fi

## ---- 5. Sanity --------------------------------------------------------------
## libR.so lives outside the default loader path; R's own launcher normally
## provides it via LD_LIBRARY_PATH, but here the binary runs straight from
## the shell.
export LD_LIBRARY_PATH="${R_HOME_DIR:-/usr/lib/R}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
## ARMPL sits outside the loader's search path; the packaged bundle carries
## it with an $ORIGIN rpath, but this in-tree binary needs to be told.
[ -n "${ARMPL_DIR:-}" ] && export LD_LIBRARY_PATH="$ARMPL_DIR/lib:$LD_LIBRARY_PATH"
"$PREFIX/bin/inla" -ping
echo "OK: inla built and installed in $PREFIX/bin"
