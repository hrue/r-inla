#!/usr/bin/env bash
# build-windows.sh — native Windows from-source build of the R-INLA `inla` engine.
#
# Distilled from .github/workflows/windows-native-build.yml so you can build the
# STATICALLY-LINKED engine on YOUR OWN Windows machine, with NO WSL and NO GitHub CI.
#
# Run it inside an Rtools45 MSYS2 UCRT64 shell (see BUILD-WINDOWS.md). It:
#   1. installs the UCRT64 toolchain + numeric-lib deps via pacman,
#   2. vendors the header-only SIMDE dependency,
#   3. builds + installs GMRFLib (and copies the vendored taucs static lib),
#   4. generates a MinGW import lib for R.dll and builds the Rmath forwarder DLL,
#   5. generates the cgeneric headers,
#   6. builds inlaprog with static numeric/runtime libs (R + Rmath stay dynamic),
#   7. assembles the drop-in bundle under ./bundle and prints the install step.
#
# The result matches what the CI produces: inla.exe + Rmathfwd.dll (+ libcrypto-3-x64.dll,
# unless you set STATIC_CRYPTO=1) — everything else baked into inla.exe.
set -euo pipefail

# --- 0. sanity: must be an MSYS2 UCRT64 shell with R on PATH ----------------------------
if [ "${MSYSTEM:-}" != "UCRT64" ]; then
  echo "ERROR: run this from the 'Rtools45 UCRT64' / 'MSYS2 UCRT64' shell (MSYSTEM=UCRT64)." >&2
  echo "       Current MSYSTEM='${MSYSTEM:-unset}'." >&2
  exit 1
fi
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not on PATH. Add your R 'bin' dir to PATH." >&2; exit 1; }

# Run from the repo root (this script lives in build-configs/windows-native/).
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
echo "== repo root: $REPO_ROOT"

STATIC_CRYPTO="${STATIC_CRYPTO:-1}"   # 1 = static libcrypto (2-file bundle, CI-proven); 0 = dynamic crypto
PREFIX="$PWD/local"

# --- 1. install UCRT64 deps (idempotent) -----------------------------------------------
echo "== installing UCRT64 dependencies via pacman ..."
pacman -S --needed --noconfirm \
  base-devel \
  mingw-w64-ucrt-x86_64-gcc \
  mingw-w64-ucrt-x86_64-gcc-fortran \
  mingw-w64-ucrt-x86_64-make \
  mingw-w64-ucrt-x86_64-pkgconf \
  mingw-w64-ucrt-x86_64-openblas \
  mingw-w64-ucrt-x86_64-gsl \
  mingw-w64-ucrt-x86_64-metis \
  mingw-w64-ucrt-x86_64-muparser \
  mingw-w64-ucrt-x86_64-zlib \
  mingw-w64-ucrt-x86_64-openssl \
  mingw-w64-ucrt-x86_64-suitesparse \
  mingw-w64-ucrt-x86_64-libtool \
  mingw-w64-ucrt-x86_64-cmake

# --- 2a. vendor SIMDE (header-only; no MSYS2 package exists) ----------------------------
if [ ! -f "$PWD/extern/simde/simde/x86/sse2.h" ]; then
  echo "== vendoring SIMDE v0.8.2 ..."
  rm -rf "$PWD/extern/simde"
  git clone --depth 1 --branch v0.8.2 https://github.com/simd-everywhere/simde.git "$PWD/extern/simde"
fi
SIMDE_INC="-I$PWD/extern/simde"

# --- 2b. vendor + build STATIC muParser (UCRT64 ships ONLY the shared lib) --------------
# Without this, -lmuparser has no static .a; a dynamic libmuparser.dll would drag
# libstdc++-6.dll/libgcc_s/libwinpthread back into the bundle. Build our own static lib.
if [ ! -d "$PWD/extern/muparser" ]; then
  echo "== vendoring + building static muParser v2.3.5 ..."
  git clone --depth 1 --branch v2.3.5 https://github.com/beltoforion/muparser.git "$PWD/extern/muparser"
fi
# CMAKE_POLICY_VERSION_MINIMUM=3.5 overrides muParser's old cmake_minimum_required,
# which CMake 4.x otherwise rejects.
cmake -S "$PWD/extern/muparser" -B "$PWD/extern/muparser/build" \
  -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DBUILD_SHARED_LIBS=OFF -DENABLE_SAMPLES=OFF -DENABLE_OPENMP=OFF
cmake --build "$PWD/extern/muparser/build" -j
MUPARSER_A="$(find "$PWD/extern/muparser/build" -name 'libmuparser*.a' | head -1)"
test -n "$MUPARSER_A" && test -f "$MUPARSER_A" || { echo "ERROR: static libmuparser.a not built" >&2; exit 1; }
MUPARSER_INC="-I$PWD/extern/muparser/include"
echo "== static muParser: $MUPARSER_A"

# --- 3. build + install GMRFLib, stage taucs -------------------------------------------
echo "== building GMRFLib ..."
mkdir -p "$PREFIX/lib" "$PREFIX/include" "$PREFIX/bin"
make -C gmrflib \
  PREFIX="$PREFIX" LEXTPREFIX=/ucrt64 \
  CC=gcc CXX=g++ FC=gfortran \
  FLAGS="-O2 -fopenmp -pthread -DINLA_WITH_OPENBLAS -DINLA_WITH_SIMDE $SIMDE_INC"
make -C gmrflib PREFIX="$PREFIX" LEXTPREFIX=/ucrt64 install
cp -v gmrflib/taucs/libtaucs.a "$PREFIX/lib/" || true

# --- 4. R import lib + Rmath forwarder DLL ---------------------------------------------
echo "== preparing R import lib + Rmath forwarder DLL ..."
RHOME="$(cygpath -u "$(Rscript -e 'cat(normalizePath(R.home()))' | tr -d '\r')")"
RBIN="$RHOME/bin/x64"; [ -d "$RBIN" ] || RBIN="$RHOME/bin"
echo "   RHOME=$RHOME"
echo "   RBIN=$RBIN"
(
  cd "$RBIN"
  gendef R.dll >/dev/null 2>&1
  dlltool -d R.def -l libR.dll.a && echo "   generated libR.dll.a"
  # INLA is built MATHLIB_STANDALONE and references Rmath functions by their UNPREFIXED
  # names (qgamma, dnorm4, ...). Windows R ships no libRmath, but R.dll exports the same
  # functions Rf_-prefixed. Build a forwarder DLL: unprefixed export `qgamma` FORWARDS at
  # load time to R.Rf_qgamma. (A plain import-lib alias does NOT work.)
  echo "LIBRARY Rmathfwd.dll" > Rmathfwd.def
  echo "EXPORTS" >> Rmathfwd.def
  grep -oE 'Rf_[A-Za-z0-9_]+' R.def | sort -u > rf-exports.txt
  grep -oE '\b(d|p|q|r)(norm|unif|gamma|beta|binom|pois|t|f|chisq|cauchy|exp|geom|hyper|lnorm|logis|nbinom|weibull|wilcox|signrank|nbeta|nchisq|nt|nf|nbinom_mu)[0-9]*\b|\b(gammafn|lgammafn|digamma|trigamma|tetragamma|pentagamma|beta|lbeta|choose|lchoose|bessel_i|bessel_j|bessel_k|bessel_y|psigamma|sign|fsign|fround|fprec|ftrunc|log1p|expm1|bessel_i_ex|bessel_k_ex)\b' "$RHOME/include/Rmath.h" | sort -u > rmath-names.txt
  n=0
  while read fn; do
    if grep -qx "Rf_$fn" rf-exports.txt; then echo "$fn = R.Rf_$fn" >> Rmathfwd.def; n=$((n+1)); fi
  done < rmath-names.txt
  echo "   Rmath forwarder: $n functions forwarded to R.dll Rf_ exports"
  printf 'int __rmathfwd_dummy;\n' > rmathfwd_dummy.c
  gcc -c rmathfwd_dummy.c -o rmathfwd_dummy.o
  gcc -shared -o Rmathfwd.dll rmathfwd_dummy.o Rmathfwd.def -Wl,--out-implib,libRmathfwd.dll.a
  ls -la Rmathfwd.dll libRmathfwd.dll.a
)

# --- 5. cgeneric headers ---------------------------------------------------------------
echo "== generating cgeneric headers ..."
( cd external-packages && bash ./update-cgeneric )

# --- 6. build inlaprog, STATICALLY linked ----------------------------------------------
# -static           : prefer libfoo.a over libfoo.dll.a for every -l, incl. the compiler
#                     runtime the driver appends (libgomp/libwinpthread/libstdc++/libgcc).
#                     Libs that only ship an import lib (R, Rmathfwd) stay dynamic.
# RLIB_LIB          : force -lRmathfwd/-lR dynamic even under -static (-Wl,-Bdynamic island).
# EXTLIBS2          : numeric libs static; crypto dynamic by default (openssl static drags
#                     in ws2_32/crypt32); set STATIC_CRYPTO=1 to try static (adds those libs).
# EXTLIBS3          : Fortran runtime static; -lm is a system stub.
echo "== building inlaprog (static link) ; STATIC_CRYPTO=$STATIC_CRYPTO"
# muParser is linked as our full-path static archive ($MUPARSER_A). crypto:
#   STATIC_CRYPTO=1 -> static libcrypto + its Windows import libs (default; smallest bundle)
#   STATIC_CRYPTO=0 -> dynamic libcrypto (ships libcrypto-3-x64.dll) via a -Bdynamic island
if [ "$STATIC_CRYPTO" = "1" ]; then
  CRYPTO="-lcrypto -lws2_32 -lcrypt32 -lbcrypt -ladvapi32 -luser32"
else
  CRYPTO="-Wl,-Bdynamic -lcrypto -Wl,-Bstatic"
fi
make -C inlaprog \
  PREFIX="$PREFIX" LEXTPREFIX=/ucrt64 \
  CC=gcc CXX=g++ FC=gfortran \
  FLAGS="-std=gnu99 -O2 -fopenmp -pipe -DINLA_WITH_OPENBLAS -DINLA_WITH_SIMDE -DINLA_WITH_MUPARSER -DMUPARSER_STATIC $SIMDE_INC $MUPARSER_INC" \
  LDFLAGS="-O2 -fopenmp -pipe -static -static-libgcc -static-libstdc++" \
  RLIB_INC="-DINLA_WITH_LIBR -I$RHOME/include" \
  RLIB_LIB="-L$RBIN -Wl,-Bdynamic -lRmathfwd -lR -Wl,-Bstatic" \
  EXTLIBS1="-L$PREFIX/lib -lGMRFLib -ltaucs" \
  EXTLIBS2="-lgsl -lmetis -lopenblas $MUPARSER_A -lz -lltdl $CRYPTO" \
  EXTLIBS3="-lgfortran -lquadmath -lm"

test -f inlaprog/inla.exe || { echo "ERROR: inla.exe was not produced" >&2; exit 1; }
echo "== built inlaprog/inla.exe"
echo "== ldd (dependency check) — only R.dll/Rmathfwd.dll (+ maybe libcrypto) should be non-system:"
ldd inlaprog/inla.exe | sort || true

# --- 7. assemble drop-in bundle --------------------------------------------------------
echo "== assembling ./bundle ..."
rm -rf bundle && mkdir -p bundle
cp inlaprog/inla.exe bundle/
cp "$RBIN/Rmathfwd.dll" bundle/
deps_of() { ldd "$1" 2>/dev/null | awk '{print $1}'; }
declare -A seen
queue="$(deps_of inlaprog/inla.exe)"
while [ -n "$queue" ]; do
  next=""
  for name in $queue; do
    [ -n "${seen[$name]:-}" ] && continue
    seen[$name]=1
    src="/ucrt64/bin/$name"
    if [ -f "$src" ]; then cp -n "$src" bundle/ && echo "   bundled $name"; next="$next $(deps_of "$src")"; fi
  done
  queue="$next"
done

echo
echo "=== BUNDLE READY: $PWD/bundle ==="
ls -la bundle/
echo
echo "Install into an existing INLA (run in R):"
echo '    d <- system.file("bin/windows/64bit", package = "INLA")'
echo "    file.copy(list.files(\"$(cygpath -m "$PWD/bundle")\", full.names = TRUE), d, overwrite = TRUE)"
echo
echo "The default engine (inla.call.builtin()) then resolves to this inla.exe."
