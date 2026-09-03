#!/usr/bin/env bash
## Build inla on macOS (Apple Silicon) with the Homebrew GCC toolchain,
## following the same stages as ci/build.sh. Differences from Linux, per
## the upstream macOS recipe: ld64's -force_load instead of
## --whole-archive, no NUMA/hwloc/CLONE_TARGETS (Linux-only), BLAS/LAPACK
## from R's own libRblas/libRlapack, static libstdc++/libgcc, and rpaths
## into R's lib dir so the binary runs without environment setup.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
## Shared versions (MIMALLOC_VERSION and friends); build.sh has sourced this
## from the start, this script never did, so any key used here was silently
## empty.
[ -f "$ROOT/ci/toolchain.env" ] && . "$ROOT/ci/toolchain.env"
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
## Version reported by the built binary. It is the R package's Version from
## rinla/DESCRIPTION, so `inla -V` and packageVersion("INLA") agree for anyone
## who installs the R package and the binary from the same commit. The short
## commit is kept in INLA_TAG (a string, shown by `inla -v` as "Build tag"),
## which is where the build-traceability belongs; GITCOMMIT has to stay a bare
## preprocessor token, so it carries the version alone.
SHA=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
TAG=$(sed -n 's/^Version:[[:space:]]*//p' "$ROOT/rinla/DESCRIPTION" 2>/dev/null | head -1)
[ -n "$TAG" ] || TAG=$SHA
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
 -DGITCOMMIT=$TAG -DINLA_TAG='\"$TAG ($SHA)\"' \
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

## BLAS/LAPACK: Apple Silicon uses ARMPL, Intel a SERIAL OpenBLAS (see the
## reasoning at the x86_64 branch below; Accelerate remains only as the
## fallback when that OpenBLAS was not built). Either falls back to R's own
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
        ## libastring is left out on purpose: see ci/build.sh for the
        ## out-of-bounds memchr in ARMPL's own cpuinfo parsing.
        ## ARMPL's archive carries C++ built against libc++ (the std::__1
        ## symbols), while this binary is linked by g++ against libstdc++.
        ## Both can coexist -- they use different inline namespaces and
        ## ARMPL's C++ never crosses the boundary -- but libc++ has to be on
        ## the link line or those symbols stay undefined.
        BLASLIBS="$ARMPL_A $AMATH -L$ARMPL_DIR/lib -lc++"
    else
        echo "BLAS: ARMPL (shared) at $ARMPL_DIR"
        BLASLIBS="-L$ARMPL_DIR/lib -Wl,-rpath,$ARMPL_DIR/lib -larmpl -lamath"
    fi
elif [ "$ARCH" = x86_64 ]; then
    ## NOT Accelerate: its internal threading has no reliable control (no
    ## runtime API; VECLIB_MAXIMUM_THREADS is legacy and partially ignored),
    ## and under INLA's nested parallelism an uncontrollable BLAS pool
    ## oversubscribes -- upstream observed exactly such timing anomalies.
    ## A SERIAL locking-safe OpenBLAS (built by ci/deps-macos.sh with the
    ## same recipe as the sTiles Intel lane) never threads itself, so all
    ## parallelism stays where INLA puts it. Same policy as every other
    ## platform: sequential MKL on Linux x86_64, serial ARMPL on arm64.
    OB=${OPENBLAS_SERIAL_ROOT:-$ROOT/local/openblas}
    if [ -f "$OB/lib/libopenblas.dylib" ]; then
        echo "BLAS: OpenBLAS (serial, locking-safe) at $OB"
        BLASLIBS="-L$OB/lib -Wl,-rpath,$OB/lib -lopenblas"
    else
        echo "BLAS: Accelerate framework (serial OpenBLAS not found at $OB)"
        FLAGS="$FLAGS -DINLA_WITH_FRAMEWORK_ACCELERATE"
        BLASLIBS="-framework Accelerate"
    fi
else
    echo "BLAS: R's Rblas/Rlapack"
    BLASLIBS="-lRblas -lRlapack"
fi

## sTiles: an alternative sparse-matrix backend the sources already support
## (gmrflib/smtp-stiles.c, guarded by INLA_WITH_STILES). Only a prebuilt
## library and its header are needed -- ci/fetch-stiles.sh stages them from
## a published release, so no sTiles source is involved.
##
## -L goes AFTER the Homebrew paths on the link line below, deliberately:
## the macOS libstiles ships its own copies of the GCC runtime
## (libgomp/libgfortran/libstdc++) beside itself, and this build must keep
## resolving those from Homebrew. Nothing here should let the linker prefer
## the copies that travel with the solver.
WITH_STILES=${WITH_STILES:-0}
STILES_LIBS=""
STILES_INC=""
if [ "$WITH_STILES" = 1 ]; then
    STILES_DIR=${STILES_DIR:-$PREFIX/stiles}
    [ -f "$STILES_DIR/include/stiles.h" ] \
        || { echo "ERROR: no stiles.h under $STILES_DIR (run ci/fetch-stiles.sh)"; exit 1; }
    echo "== sTiles: $STILES_DIR =="
    ## -ltileindexer is NOT needed: released libraries carry it inside.
    STILES_INC="-DINLA_WITH_STILES -I$STILES_DIR/include"
    STILES_LIBS="-L$STILES_DIR/lib -Wl,-rpath,$STILES_DIR/lib -lstiles"
fi
FLAGS="$FLAGS $STILES_INC"

## ---- 1. External model packages --------------------------------------------
for d in "$EPATH"/*/; do
    [ -f "$d/Makefile" ] && make -C "$d" clean >/dev/null
done
rm -f "$EPATH"/lib*.a
## -I for libtool: defining INLA_WITH_EXTERNAL_PACKAGES (in INC below) is what
## makes an external package's header include <ltdl.h>. Linux picks that up
## from /usr/include with no -I at all, so the omission stayed invisible until
## a package started honouring the macro. Homebrew's libtool is keg-only, i.e.
## NOT symlinked into $BREW/include, which is why the link line further down
## already carries -L$BREW/opt/libtool/lib; the compile side needs the matching
## include path or ltdl.h is not found.
##
## Keep comments OUT of the command below: every line up to the closing paren
## is one continued command, so a "#" line there is not a comment, it is an
## argument to ./build and is forwarded into make.
( cd "$EPATH" && ./build \
      CC="$CC" CXX="$CXX" FC="$FC" \
      FLAGS="" \
      INC="-DINLA_WITH_EXTERNAL_PACKAGES -I$BREW/include/eigen3 \
           -I$BREW/opt/libtool/include \
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
               -static-libstdc++ -static-libgcc -lm -ldl $STILES_LIBS" \
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
## ---- mimalloc: shipped beside the binary, never linked (see build.sh) -------
if [ -n "${MIMALLOC_VERSION:-}" ]; then
    rm -rf "$ROOT/mimalloc-src"
    git clone -q --depth 1 --branch "$MIMALLOC_VERSION" \
        https://github.com/microsoft/mimalloc "$ROOT/mimalloc-src"
    ## MI_MALLOC_OVERRIDE is what makes the library actually REPLACE
    ## malloc/free when preloaded; cmake sets it, the single-TU build must
    ## too, or the .so exports no allocator entry points at all (verified
    ## with nm before this flag was added).
    $CC -O2 -DNDEBUG -DMI_MALLOC_OVERRIDE -fPIC -shared -pthread \
        -I"$ROOT/mimalloc-src/include" \
        "$ROOT/mimalloc-src/src/static.c" \
        -install_name @rpath/libmimalloc.dylib \
        -o "$PREFIX/lib/libmimalloc.dylib" \
        || { echo "ERROR: mimalloc $MIMALLOC_VERSION failed to build"; exit 1; }
    echo "mimalloc: $MIMALLOC_VERSION -> $PREFIX/lib/libmimalloc.dylib"
fi

BLAS_DESC=${BLAS_DESC:-"$BLASLIBS"}
bash "$ROOT/ci/write-buildinfo.sh" "$PREFIX/BUILDINFO" "$CC" "$FLAGS" "$BLAS_DESC"
"$PREFIX/bin/inla" -ping
echo "OK: inla built and installed in $PREFIX/bin (macOS arm64)"
