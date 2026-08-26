#!/usr/bin/env bash
## Fetch a prebuilt libstiles for this platform and stage it as
## <prefix>/stiles/{lib,include}, which is what ci/build.sh expects when
## WITH_STILES=1.
##
## The release assets are public, so no credentials are involved and no
## sTiles source is needed: a binary and stiles.h are the whole interface.
##
## Knobs:
##   STILES_REPO     owner/name publishing the assets
##   STILES_TAG      a release tag, or "latest" (default)
##   STILES_ASSET    override the asset name for this platform
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
DEST=$PREFIX/stiles
STILES_REPO=${STILES_REPO:-esmail-abdulfattah/sTiles}
STILES_TAG=${STILES_TAG:-latest}

## Which asset matches this machine. The GCC-built variants are the ones an
## INLA process can link: it is built with GCC, so it uses GNU libgomp, and
## a clang-built library would bring LLVM's libomp into the same process.
if [ -z "${STILES_ASSET:-}" ]; then
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)   STILES_ASSET=libstiles-linux-x86_64.zip ;;
        Linux-aarch64)  STILES_ASSET=libstiles-linux-arm64.zip ;;
        Darwin-arm64)   STILES_ASSET=libstiles-macos-arm64-gcc-armpl.zip
                        STILES_ASSET_FALLBACK=libstiles-macos-apple-arm64-armpl.zip ;;
        Darwin-x86_64)  STILES_ASSET=libstiles-macos-intel-x86_64-gcc-openblas.zip ;;
        *)              echo "ERROR: no libstiles asset defined for $(uname -s)-$(uname -m)"; exit 1 ;;
    esac
fi

if [ "$STILES_TAG" = latest ]; then
    URL="https://github.com/$STILES_REPO/releases/latest/download/$STILES_ASSET"
else
    URL="https://github.com/$STILES_REPO/releases/download/$STILES_TAG/$STILES_ASSET"
fi

echo "== fetching $STILES_ASSET from $STILES_REPO ($STILES_TAG) =="
rm -rf "$DEST" /tmp/stiles-dl
mkdir -p "$DEST/lib" "$DEST/include" /tmp/stiles-dl
if ! EFFECTIVE=$(curl -fL --retry 3 --retry-delay 10 -o /tmp/stiles-dl/asset.zip \
                 -w '%{url_effective}' "$URL"); then
    ## older releases may carry the asset under its previous name
    if [ -n "${STILES_ASSET_FALLBACK:-}" ]; then
        echo "  $STILES_ASSET not in this release; trying $STILES_ASSET_FALLBACK"
        URL=${URL/$STILES_ASSET/$STILES_ASSET_FALLBACK}
        EFFECTIVE=$(curl -fL --retry 3 --retry-delay 10 -o /tmp/stiles-dl/asset.zip \
                     -w '%{url_effective}' "$URL") \
            || { echo "ERROR: cannot download $URL"; exit 1; }
    else
        echo "ERROR: cannot download $URL"; exit 1
    fi
fi
## Record which release "latest" resolved to; the release notes table
## reads this. The effective URL usually ends on the CDN (no tag in it),
## so ask the /releases/latest redirect for the tag instead. Never fatal:
## a miss records the requested tag.
RELTAG=$(echo "$EFFECTIVE" | grep -oE '/download/[^/]+/' | head -1 | sed 's|/download/||; s|/||' || true)
if [ -z "$RELTAG" ] && [ "$STILES_TAG" = latest ]; then
    RELTAG=$(curl -sI "https://github.com/$STILES_REPO/releases/latest" 2>/dev/null \
             | grep -i '^location:' | grep -oE '/tag/[^[:space:]]+' | sed 's|/tag/||' | tr -d '\r' || true)
fi
echo "${RELTAG:-$STILES_TAG}" > "$DEST/RELEASE"
echo "  release: $(cat "$DEST/RELEASE")"
## unzip is not in every container image this runs in (the manylinux ones
## ship neither unzip nor bsdtar), so fall back to python's zipfile rather
## than adding a package to four different dependency scripts.
if command -v unzip >/dev/null 2>&1; then
    ( cd /tmp/stiles-dl && unzip -q asset.zip )
elif command -v bsdtar >/dev/null 2>&1; then
    ( cd /tmp/stiles-dl && bsdtar xf asset.zip )
elif command -v python3 >/dev/null 2>&1; then
    python3 -c 'import sys,zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])' \
            /tmp/stiles-dl/asset.zip /tmp/stiles-dl
else
    echo "ERROR: no way to unpack a zip here (unzip, bsdtar and python3 all missing)"; exit 1
fi

## The layout inside the asset has moved between releases, so locate the
## pieces rather than assume where they are.
HDR=$(find /tmp/stiles-dl -name 'stiles.h' | head -1)
LIB=$(find /tmp/stiles-dl \( -name 'libstiles.so*' -o -name 'libstiles.dylib' \
                             -o -name 'libstiles.dll' -o -name 'stiles.dll' \) | head -1)
[ -n "$HDR" ] || { echo "ERROR: stiles.h not in the asset"; find /tmp/stiles-dl | head -20; exit 1; }
[ -n "$LIB" ] || { echo "ERROR: no libstiles in the asset"; find /tmp/stiles-dl | head -20; exit 1; }
cp "$HDR" "$DEST/include/"
## Everything beside the library travels too: a portable libstiles may ship
## its own dependencies next to it.
cp -a "$(dirname "$LIB")"/. "$DEST/lib/"

echo "  header:  $DEST/include/stiles.h"
echo "  library: $(ls "$DEST"/lib/libstiles.* 2>/dev/null | head -1)"
ls -lh "$DEST/lib" | head -10

## ---------------------------------------------------------------------------
## Same compiler on both sides, checked rather than assumed.
##
## The shared toolchain file pins the compiler where a package manager lets us
## (GCC_PREFER for the modern Linux lanes, GCC_TOOLSET inside the manylinux
## containers). It cannot on macOS or Windows: Homebrew and MSYS2 serve only
## their current package, so two lanes built on different days get different
## point releases. Pinning is therefore not enough, and the pairing has already
## drifted twice without anything saying so: the baseline Linux bundle carried a
## library built by GCC 14 inside a binary built by GCC 15, and the Windows DLL
## was built by MinGW 16.2.0 while inla.exe used 16.1.1.
##
## So compare what actually built each side. The asset records its compiler in
## BUILDINFO; this build knows its own. A mismatch fails here, where it is one
## line to read, instead of surviving into a release nobody questions.
## STILES_ALLOW_COMPILER_MISMATCH=1 downgrades it to a warning for the case
## where a deliberate mismatch is being tested.
BI=$(find /tmp/stiles-dl -maxdepth 3 -name BUILDINFO | head -1)
if [ -n "$BI" ]; then
    ## "g++ (GCC) 14.2.1 ..." / "gcc-16 (Ubuntu ...) 16.0.1 ..." / MSYS2 rev
    lib_cc=$(grep -i '^compiler' "$BI" | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    ## Compare against the compiler that builds for the SAME target. The
    ## Windows lane cross-builds from Linux, so $CXX there is the host g++ and
    ## comparing it with a DLL built by MSYS2 would fail for no reason.
    ## Resolve the compiler the way ci/build.sh does (newest versioned g++,
    ## else the default): the Ubuntu PPA installs g++-16 WITHOUT repointing
    ## bare g++, so querying g++ here compared the library against a compiler
    ## the build never uses -- 16.0.1 vs 13.3.0 on a lane that builds with 16.
    if [ -z "${CXX:-}" ]; then
        ## manylinux first: the pinned gcc-toolset lives under /opt/rh, not as
        ## a versioned /usr/bin binary, and the image leaves the OLD toolset on
        ## PATH at fetch time -- bare g++ answered 14.2.1 on a lane whose build
        ## uses toolset-15.
        HERE=$(cd "$(dirname "$0")" && pwd)
        [ -f "$HERE/toolchain.env" ] && . "$HERE/toolchain.env"
        if [ -n "${GCC_TOOLSET:-}" ] && [ -x "/opt/rh/gcc-toolset-$GCC_TOOLSET/root/usr/bin/g++" ]; then
            CXX="/opt/rh/gcc-toolset-$GCC_TOOLSET/root/usr/bin/g++"
        elif [ -n "${GCC_TOOLSET:-}${GCC_PREFER:-}" ]; then
            ## A toolchain IS pinned (by GCC_TOOLSET or GCC_PREFER) but
            ## nothing upstream of this point resolved it to an actual
            ## compiler: no exported CXX, no matching /opt/rh toolset.
            ## Globbing /usr/bin for "whatever's highest" here is exactly
            ## how this checker itself would end up comparing the library
            ## against a compiler the build never uses. Fail loudly instead;
            ## fix the step that was supposed to install and export the pin.
            echo "ERROR: a compiler is pinned (GCC_TOOLSET=${GCC_TOOLSET:-} GCC_PREFER=${GCC_PREFER:-})" >&2
            echo "       but this step has no CXX and no matching gcc-toolset." >&2
            echo "       The step that installs it should export CC/CXX/FC" >&2
            echo "       (see ci/deps-ubuntu.sh) or run before this one." >&2
            exit 1
        else
            ## genuinely unpinned (canary / newest-gcc lanes): best effort
            CXX=$(ls /usr/bin/g++-1[0-9] 2>/dev/null | sort -V | tail -1 || true)
        fi
    fi
    own_cxx=${CXX:-g++}
    if [ -f "$DEST/lib/libstiles.dll" ]; then
        for c in ${MINGW_CXX:-} x86_64-w64-mingw32ucrt-g++ x86_64-w64-mingw32-g++; do
            command -v "$c" >/dev/null 2>&1 && { own_cxx=$c; break; }
        done
    fi
    own_cc=$($own_cxx --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    echo "  compiler: library $lib_cc, this build ${own_cc:-unknown} ($own_cxx)"
    if [ -n "$lib_cc" ] && [ -n "$own_cc" ] && [ "$lib_cc" != "$own_cc" ]; then
        if [ "${STILES_ALLOW_COMPILER_MISMATCH:-0}" = "1" ]; then
            echo "  WARNING: compiler mismatch ($lib_cc vs $own_cc), allowed by request"
        else
            echo "ERROR: libstiles was built by GCC $lib_cc but this build uses $own_cc." >&2
            echo "       The bundle would ship a library compiled by a different" >&2
            echo "       compiler than the binary around it. Align them through the" >&2
            echo "       shared toolchain file (GCC_PREFER / GCC_TOOLSET), or rebuild" >&2
            echo "       the sTiles asset. Set STILES_ALLOW_COMPILER_MISMATCH=1 to" >&2
            echo "       proceed anyway." >&2
            exit 1
        fi
    fi
fi

## The library is portable but not self-contained: it expects a few system
## libraries (libnuma, libhwloc, libgfortran) from the host. Resolve it now,
## so a missing one is a clear message here rather than a link or load
## failure later. Its OpenMP runtime is reported for the same reason: an
## INLA process is GCC-built and carries libgomp, and a second runtime in
## the same process is the classic "OMP: Error #15" abort.
## A Windows asset is fetched from a Linux cross-build host, so nothing
## here can load it; the bare-bundle test on a real Windows runner is what
## checks that one.
if [ -f "$DEST/lib/libstiles.dll" ]; then
    echo "== Windows DLL (cross-build host cannot inspect it) =="
    ls -1 "$DEST/lib" | sed 's/^/  /'
    exit 0
fi

case "$(uname -s)" in
    Linux)
        SO=$(ls "$DEST"/lib/libstiles.so* 2>/dev/null | head -1)
        echo "== runtime dependencies =="
        ldd "$SO" | sed 's/^/  /'
        if ldd "$SO" | grep -q "not found"; then
            echo "ERROR: libstiles needs libraries this machine does not have (see above)"
            exit 1
        fi
        echo "  OpenMP runtime: $(ldd "$SO" | grep -oE 'lib(gomp|iomp5|omp)\.so[^ ]*' | head -1 || echo none)"
        ;;
    Darwin)
        DY=$(ls "$DEST"/lib/libstiles.dylib 2>/dev/null | head -1)
        echo "== runtime dependencies =="
        otool -L "$DY" | sed 's/^/  /'
        ;;
esac
