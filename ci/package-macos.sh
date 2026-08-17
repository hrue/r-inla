#!/usr/bin/env bash
## Bundle the macOS binary with every non-system dylib it needs, so it runs
## on a Mac that has neither Homebrew nor ARMPL installed. Apple's own
## frameworks (Accelerate, libSystem) stay external: they are part of the
## operating system.
##
## R is the deliberate exception when the binary was linked against libR
## (WITH_LIBR=1): it comes from the user's own R installation, exactly as
## upstream's macOS binaries expect. A WITH_LIBR=2 build has no libR
## dependency at all.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-$ROOT/local}
OUT=$ROOT/dist
BIN=$PREFIX/bin/inla
ARCH=$(uname -m)

rm -rf "$OUT"
mkdir -p "$OUT/bin" "$OUT/lib"
cp "$BIN" "$OUT/bin/inla"

echo "== dependencies before bundling =="
otool -L "$OUT/bin/inla" || true

## dylibbundler copies every non-system dependency (recursively) next to
## the binary and rewrites the load commands to @loader_path.
## </dev/null: it prompts interactively when it cannot resolve something,
## which would otherwise hang the job.
dylibbundler --overwrite-files --bundle-deps --create-dir \
    --fix-file "$OUT/bin/inla" \
    --dest-dir "$OUT/lib" \
    --install-path @loader_path/../lib/ </dev/null

## dylibbundler can leave a duplicate @loader_path rpath, which dyld on
## recent macOS rejects outright; the load commands already carry the full
## @loader_path/../lib/ path, so no rpath entry is needed at all. Every
## edited file has to be re-signed (ad-hoc) on Apple Silicon.
for f in "$OUT/bin/inla" "$OUT"/lib/*.dylib; do
    [ -f "$f" ] || continue
    while install_name_tool -delete_rpath @loader_path/ "$f" 2>/dev/null; do :; done
    codesign -f -s - "$f" >/dev/null 2>&1 || true
done

echo "== dependencies after bundling =="
otool -L "$OUT/bin/inla"
ls -lh "$OUT/lib" 2>/dev/null || true

## Nothing outside the system (and the user's R) may remain referenced by
## absolute path, otherwise the artifact only runs on this build machine.
if otool -L "$OUT/bin/inla" "$OUT"/lib/*.dylib 2>/dev/null \
     | grep -E '/opt/homebrew|/usr/local/(opt|Cellar)|/opt/arm|/opt/local'; then
    echo "ERROR: a third-party library is still referenced by absolute path"
    exit 1
fi

echo "== BLAS in the bundle =="
if ls "$OUT"/lib 2>/dev/null | grep -iE 'armpl|amath|astring|blas'; then
    echo "  ^ shipped beside the binary"
else
    n=$(nm -U "$OUT/bin/inla" 2>/dev/null | grep -ciE 'dgemm|dpotrf' || true)
    echo "  embedded in the binary, or Accelerate (BLAS symbols defined: $n)"
fi

"$OUT/bin/inla" -ping

NAME=inla-macos-$ARCH-portable
tar -C "$OUT" -czf "$ROOT/$NAME.tar.gz" bin lib
echo "OK: portable bundle $NAME.tar.gz ($(du -h "$ROOT/$NAME.tar.gz" | cut -f1))"
