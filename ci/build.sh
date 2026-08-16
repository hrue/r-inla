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
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
WITH_LIBR=${WITH_LIBR:-1}
JOBS=${JOBS:-$(nproc)}
CC=${CC:-gcc}
CXX=${CXX:-g++}
FC=${FC:-gfortran}
EPATH=$ROOT/external-packages

## Version string compiled into the binary (shown by inla -V and inla -ping);
## falls back to "devel" outside a git checkout.
TAG=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo devel)

echo "== building $TAG with CC=$CC CXX=$CXX FC=$FC ($($CC --version | head -1)) =="

## The devel feature set, minus PARDISO (requires a license; libpardiso.c
## provides stubs so the link still closes) and minus MKL (this build links
## OpenBLAS). INLA_WITH_SIMDE needs the SIMDE headers (libsimde-dev).
## No -march: baseline x86-64, with the fast per-CPU paths supplied at run
## time by the INLA_WITH_CLONE_TARGETS function clones.
FLAGS="-O2 -mtune=generic -ftree-vectorize -funroll-loops -pipe -pthread \
 -fopenmp -fopenmp-simd -flax-vector-conversions \
 -DINLA_WITH_SIMDE -DINLA_WITH_DEVEL -DINLA_WITH_CLONE_TARGETS \
 -DINLA_WITH_EXTERNAL_PACKAGES -DINLA_WITH_MUPARSER -DINLA_WITH_NUMA \
 -DGITCOMMIT=$TAG"

## R linkage: WITH_LIBR=1 embeds the R interpreter (rgeneric models work;
## needs the shared libR that Debian/Ubuntu R provides). WITH_LIBR=0 links
## only the standalone Rmath library; everything but rgeneric works.
if [ "$WITH_LIBR" = 1 ]; then
    ## libR's directory differs per family (Debian: /usr/lib/R/lib,
    ## Fedora: /usr/lib64/R/lib); ask R itself.
    R_HOME_DIR=${R_HOME_DIR:-$(R RHOME 2>/dev/null || echo /usr/lib/R)}
    RLIB_INC="-DINLA_WITH_LIBR -I/usr/include/R -I/usr/share/R/include"
    RLIB_LIB="-lRmath -L$R_HOME_DIR/lib -lR"
else
    RLIB_INC="-I/usr/share/R/include"
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
      INC="-DINLA_WITH_EXTERNAL_PACKAGES -I/usr/include/eigen3 \
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
## -I$EPATH lets inla.c find the generated cgeneric-table.h.
## One space-separated line (ls would emit newlines, and a newline inside a
## make command-line variable truncates the link recipe at that point,
## leaving --whole-archive unterminated).
EXTOBJ=$(echo "$EPATH"/lib*.a)
make -C "$ROOT/inlaprog" -j"$JOBS" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" \
     EXTLIBS2="-Wl,--whole-archive $EXTOBJ -Wl,--no-whole-archive \
               -lgsl -lopenblas -lmuparser -lz -lmetis \
               -lnuma -lhwloc -lltdl -lcrypto -lgfortran -lquadmath -lm -ldl"
make -C "$ROOT/inlaprog" PREFIX="$PREFIX" \
     CC="$CC" CXX="$CXX" FC="$FC" \
     FLAGS="$FLAGS -I$EPATH" \
     RLIB_INC="$RLIB_INC" RLIB_LIB="$RLIB_LIB" install

## ---- 4. Sanity --------------------------------------------------------------
## libR.so lives outside the default loader path; R's own launcher normally
## provides it via LD_LIBRARY_PATH, but here the binary runs straight from
## the shell.
export LD_LIBRARY_PATH="${R_HOME_DIR:-/usr/lib/R}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$PREFIX/bin/inla" -ping
echo "OK: inla built and installed in $PREFIX/bin"
