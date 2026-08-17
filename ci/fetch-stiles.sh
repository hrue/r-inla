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
        Darwin-arm64)   STILES_ASSET=libstiles-macos-apple-arm64-armpl.zip ;;
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
curl -fL --retry 3 --retry-delay 10 -o /tmp/stiles-dl/asset.zip "$URL" \
    || { echo "ERROR: cannot download $URL"; exit 1; }
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
