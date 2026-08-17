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

## The portability floor: no glibc symbol newer than the container's 2.28.
floor=$(objdump -T "$OUT/bin/inla" | grep -oE 'GLIBC_2\.[0-9]+' | sort -uV | tail -1)
echo "glibc floor: ${floor:-none}"
case "${floor:-GLIBC_2.0}" in
    GLIBC_2.2[0-8]|GLIBC_2.1?|GLIBC_2.[0-9]) : ;;
    *) echo "ERROR: glibc floor $floor exceeds 2.28 (a host library leaked in)"; exit 1 ;;
esac

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
NAME=inla-linux-$ARCH${SUFFIX:+-$SUFFIX}-portable
tar -C "$OUT" -czf "$ROOT/$NAME.tar.gz" bin lib
echo "OK: portable bundle $NAME.tar.gz ($(du -h "$ROOT/$NAME.tar.gz" | cut -f1))"
