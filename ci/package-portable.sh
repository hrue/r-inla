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

rm -rf "$OUT"
mkdir -p "$OUT/bin" "$OUT/lib"
cp "$BIN" "$OUT/bin/inla"

## Everything ldd resolves, except the glibc family and the loader.
ldd "$BIN" | awk '/=>/ {print $3}' | while read -r so; do
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

## Prove the bundle runs from its own libraries.
"$OUT/bin/inla" -ping

tar -C "$OUT" -czf "$ROOT/inla-linux-x86_64-portable.tar.gz" bin lib
echo "OK: portable bundle $(du -h "$ROOT/inla-linux-x86_64-portable.tar.gz" | cut -f1)"
