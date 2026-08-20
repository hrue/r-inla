#!/usr/bin/env bash
## Bundle the freshly built inla binary with every non-glibc shared library
## it uses, producing one tarball that runs on any x86-64 Linux whose glibc
## is at least the build container's (manylinux_2_28 -> 2.28, i.e. distros
## from roughly 2018 on). glibc itself stays dynamic on purpose: that is
## what makes the floor mechanism work.
##
## The binary is built WITH_LIBR=2 (dlopen mode): it starts with any R or
## none, and rgeneric works by loading the running machine's own libR the
## first time an rgeneric model is used.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
OUT=$ROOT/dist
BIN=$PREFIX/bin/inla

## Libraries that live outside the loader's search path (ARMPL) must be
## resolvable here, otherwise ldd cannot report them and they would be
## missing from the bundle.
ARMPL_DIR=${ARMPL_DIR:-$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)}
[ -n "$ARMPL_DIR" ] && export LD_LIBRARY_PATH="$ARMPL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

rm -rf "$OUT"
mkdir -p "$OUT/bin" "$OUT/lib"
cp "$BIN" "$OUT/bin/inla"

## inla.run beside the binary: the upstream entry point is the SCRIPT (it
## preloads the bundled allocator when it finds one, then execs bin/inla).
## Canonical copy from utils/scripts; the per-platform R-package copy is the
## fallback for older checkouts.
RUN_SRC=""
for c in "$ROOT/utils/scripts/inla.run" "$ROOT/rinla/inst/bin/linux/64bit/inla.run"; do
    [ -f "$c" ] && { RUN_SRC=$c; break; }
done
[ -n "$RUN_SRC" ] || { echo "ERROR: no inla.run found to ship"; exit 1; }
cp "$RUN_SRC" "$OUT/bin/inla.run"
chmod 0755 "$OUT/bin/inla.run"

## The allocator itself: copied, never linked (see ci/build.sh).
if [ -f "$PREFIX/lib/libmimalloc.so" ]; then
    cp "$PREFIX/lib/libmimalloc.so" "$OUT/lib/"
else
    echo "ERROR: libmimalloc.so was not built (MIMALLOC_VERSION unset?)"; exit 1
fi

## Everything ldd resolves, except the glibc family and the loader.
ldd "$BIN" | awk '/=>/ {print $3}' | grep -v '^$' | while read -r so; do
    [ -f "$so" ] || continue
    case "$(basename "$so")" in
        libc.so*|libm.so*|libpthread.so*|libdl.so*|librt.so*|libresolv.so*|libmvec.so*|ld-linux*) ;;
        *) cp -v "$so" "$OUT/lib/" ;;
    esac
done

## The binary finds its bundled libraries relative to itself, wherever the
## bundle is unpacked; the libraries find each other the same way.
patchelf --set-rpath '$ORIGIN/../lib' "$OUT/bin/inla"
for so in "$OUT"/lib/*.so*; do
    patchelf --set-rpath '$ORIGIN' "$so"
done
## the allocator is dlopen/preload-only; give it the same self-relative rpath
patchelf --set-rpath '$ORIGIN' "$OUT/lib/libmimalloc.so" 2>/dev/null || true

## The portability floor: no glibc symbol newer than GLIBC_MAX (default
## 2.28, the manylinux container's; the modern gcc-16 lane builds on
## Ubuntu 24.04 and passes 2.39).
GLIBC_MAX=${GLIBC_MAX:-2.28}
floor=$(objdump -T "$OUT/bin/inla" | grep -oE 'GLIBC_2\.[0-9]+' | sort -uV | tail -1)
echo "glibc floor: ${floor:-none} (allowed: <= $GLIBC_MAX)"
worst=$(printf 'GLIBC_%s\n%s\n' "$GLIBC_MAX" "${floor:-GLIBC_2.0}" | sort -V | tail -1)
[ "$worst" = "GLIBC_$GLIBC_MAX" ] \
    || { echo "ERROR: glibc floor $floor exceeds $GLIBC_MAX (a host library leaked in)"; exit 1; }

## One OpenMP runtime in the shipped bundle as well: the check in the build
## script covers the freshly linked binary, this one covers what actually
## travels to the user.
OMPRT=$(ldd "$OUT/bin/inla" 2>/dev/null | grep -oE 'lib(gomp|iomp5|omp)\.so[^ ]*' | sort -u || true)
echo "== OpenMP runtime in the bundle: ${OMPRT:-none (static)} =="
if echo "$OMPRT" | grep -q 'libiomp5' || [ "$(echo "$OMPRT" | grep -c .)" -gt 1 ]; then
    echo "ERROR: the bundle carries an unexpected set of OpenMP runtimes"
    echo "$OMPRT"
    exit 1
fi

## Prove the bundle runs from its own libraries.
"$OUT/bin/inla" -ping
## the SHIPPED entry point must work too
"$OUT/bin/inla.run" -ping
## and the wrapper must SEE the allocator it exists to preload; both live
## inside the bundle now, so absence is a packaging bug, not a layout note.
find "$OUT/lib" -name 'libmimalloc*' | grep -q . \
    || { echo "ERROR: inla.run cannot find libmimalloc.so in $OUT/lib"; exit 1; }
echo "inla.run will preload: $(find "$OUT/lib" -name 'libmimalloc*' | tail -1)"

## The BLAS is linked statically, so it should be INSIDE the binary rather
## than beside it. Report either way: this is the one property of the
## artifact that is invisible from the outside and easy to get wrong.
echo "== BLAS in the bundle =="
if ls "$OUT"/lib | grep -iE 'blas|armpl|amath|astring'; then
    echo "  ^ shipped as shared libraries beside the binary"
else
    ## Confirm the symbols really are in the executable.
    n=$(nm -D --defined-only "$OUT/bin/inla" 2>/dev/null | grep -ciE 'dgemm|dpotrf' || true)
    echo "  embedded in the binary (BLAS symbols defined: $n)"
fi

ARCH=$(uname -m)
## ARCH_NAME overrides the uname spelling (aarch64 -> arm64) in artifact names
NAME=inla-linux-${ARCH_NAME:-$ARCH}${SUFFIX:+-$SUFFIX}-portable
cp "$PREFIX/BUILDINFO" "$OUT/BUILDINFO" 2>/dev/null || true
tar -C "$OUT" -czf "$ROOT/$NAME.tar.gz" bin lib BUILDINFO
echo "OK: portable bundle $NAME.tar.gz ($(du -h "$ROOT/$NAME.tar.gz" | cut -f1))"
