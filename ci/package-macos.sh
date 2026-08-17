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

## sTiles first, for the same reason ARMPL is handled below: its dylib is
## installed as @rpath/libstiles.dylib, and dylibbundler does not read the
## binary's rpaths, so it cannot find it.
##
## Everything beside it travels too. The macOS libstiles refers to its own
## dependencies as @loader_path/<name>, so those names must resolve in the
## same directory it ends up in. Landing them all in $OUT/lib means any
## library both sides need -- libgomp above all -- exists as ONE file there,
## which is what keeps a single OpenMP runtime in the process. (Copied
## first, so INLA's own Homebrew copies overwrite them: those are the ones
## the binary was linked against.)
STILES_DIR=${STILES_DIR:-$PREFIX/stiles}
if [ -f "$STILES_DIR/lib/libstiles.dylib" ]; then
    echo "== staging sTiles from $STILES_DIR =="
    cp -a "$STILES_DIR"/lib/*.dylib "$OUT/lib/" 2>/dev/null || true
    chmod u+w "$OUT"/lib/*.dylib
    for f in "$OUT"/lib/*.dylib; do
        install_name_tool -id "@loader_path/../lib/$(basename "$f")" "$f" 2>/dev/null || true
    done
    ## libstiles' own @loader_path/<name> references stay valid ($OUT/lib is
    ## where its siblings now live). Its @rpath ones, and the binary's
    ## @rpath/libstiles.dylib, are redirected here rather than in the ARMPL
    ## pass below: that pass only runs when ARMPL is installed, so on Intel
    ## it would never happen and the bundle would keep a reference no
    ## loader can resolve. Only names actually present in $OUT/lib are
    ## touched, leaving ARMPL's own @rpath entries for that pass.
    for f in "$OUT/bin/inla" "$OUT"/lib/*.dylib; do
        [ -f "$f" ] || continue
        for dep in $(otool -L "$f" | awk '/@rpath\// {print $1}'); do
            base=${dep#@rpath/}
            [ -f "$OUT/lib/$base" ] || continue
            install_name_tool -change "$dep" "@loader_path/../lib/$base" "$f" 2>/dev/null || true
        done
    done
    ls -1 "$OUT/lib" | sed 's/^/  /'
fi

## ARMPL next, by hand. Its libraries are referenced as @rpath/... and
## dylibbundler does not read the binary's rpaths, so it cannot find them:
## it drops into an interactive prompt and, with no input, loops on it
## forever. Copying them and rewriting the load commands ourselves removes
## the question before dylibbundler is ever run.
ARMPL_DIR=$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)
if [ -n "$ARMPL_DIR" ]; then
    echo "== resolving ARMPL from $ARMPL_DIR =="
    pending=$(otool -L "$OUT/bin/inla" | awk '/@rpath\// {print $1}' | sed 's|@rpath/||')
    while [ -n "$pending" ]; do
        next=""
        for lib in $pending; do
            [ -f "$OUT/lib/$lib" ] && continue
            src="$ARMPL_DIR/lib/$lib"
            [ -f "$src" ] || { echo "ERROR: $lib not found in $ARMPL_DIR/lib"; exit 1; }
            cp -v "$src" "$OUT/lib/"
            chmod u+w "$OUT/lib/$lib"
            ## the copy must refer to itself and to its siblings by a path
            ## relative to the binary, not through an rpath
            install_name_tool -id "@loader_path/../lib/$lib" "$OUT/lib/$lib"
            next="$next $(otool -L "$OUT/lib/$lib" | awk '/@rpath\// {print $1}' | sed 's|@rpath/||')"
        done
        pending=$(echo "$next" | tr ' ' '\n' | sort -u | tr '\n' ' ')
        [ -z "$(echo "$pending" | tr -d ' ')" ] && break
    done
    ## point every @rpath reference, in the binary and in the copies, at the
    ## bundled files
    for f in "$OUT/bin/inla" "$OUT"/lib/*.dylib; do
        [ -f "$f" ] || continue
        for dep in $(otool -L "$f" | awk '/@rpath\// {print $1}'); do
            base=${dep#@rpath/}
            install_name_tool -change "$dep" "@loader_path/../lib/$base" "$f" 2>/dev/null || true
        done
    done
fi

## Everything else (Homebrew's gsl, muparser, ltdl, crypto, gfortran...) is
## referenced by absolute path, which dylibbundler handles on its own.
## timeout + </dev/null: it must never be able to wait for input.
## macOS ships no `timeout` (it is GNU coreutils); use it only if present.
TO=$(command -v gtimeout || command -v timeout || true)
${TO:+$TO 300} dylibbundler --overwrite-files --bundle-deps --create-dir \
    --fix-file "$OUT/bin/inla" \
    --dest-dir "$OUT/lib" \
    --install-path @loader_path/../lib/ </dev/null \
  || { echo "ERROR: dylibbundler could not resolve every dependency"; exit 1; }

## Normalize any remaining ABSOLUTE third-party reference to the bundled
## copy. dylibbundler only rewrites files it discovered from the binary's
## dependency chain; the sTiles-staged dylibs arrive with their own
## internal references (the GCC runtime trio -- libgfortran, libquadmath,
## libgcc_s -- referenced as /usr/local/opt/gcc/... from the build
## machine), which nothing above touches. Every such reference whose
## basename exists in $OUT/lib is pointed there; anything else survives to
## fail the gate below, which is what an unbundlable dependency should do.
for f in "$OUT/bin/inla" "$OUT"/lib/*.dylib; do
    [ -f "$f" ] || continue
    otool -L "$f" | awk 'NR>1 && $1 ~ /^\/(usr\/local|opt\/(homebrew|arm|local))\// {print $1}'     | while read -r dep; do
        base=$(basename "$dep")
        [ -f "$OUT/lib/$base" ] || continue
        chmod u+w "$f" 2>/dev/null || true
        install_name_tool -change "$dep" "@loader_path/../lib/$base" "$f" 2>/dev/null || true
    done
done

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

## The macOS version a user actually needs is the MAXIMUM over the binary
## and every bundled library, not the deployment target we asked for: one
## dependency built for a newer OS raises the floor for the whole bundle.
## ARMPL's macOS build targets 13.0, so the Apple Silicon bundle lands
## there while the Intel one (Accelerate is a system framework) stays at
## the 11.0 we set. Printed per file so a regression is visible here rather
## than discovered by someone on an older Mac.
echo "== macOS deployment floor (ours: $MACOSX_DEPLOYMENT_TARGET) =="
FLOOR=0
for f in "$OUT/bin/inla" "$OUT"/lib/*.dylib; do
    [ -f "$f" ] || continue
    v=$(otool -l "$f" | awk '/LC_BUILD_VERSION/{f=1} f&&/minos/{print $2; exit}')
    printf '  %-38s minos %s\n' "$(basename "$f")" "${v:-?}"
    case "$v" in
        [0-9]*) [ "$(printf '%s\n%s\n' "$v" "$FLOOR" | sort -V | tail -1)" = "$v" ] && FLOOR=$v ;;
    esac
done
echo "  bundle requires macOS >= ${FLOOR:-unknown}"

"$OUT/bin/inla" -ping

NAME=inla-macos-$ARCH-portable
tar -C "$OUT" -czf "$ROOT/$NAME.tar.gz" bin lib
echo "OK: portable bundle $NAME.tar.gz ($(du -h "$ROOT/$NAME.tar.gz" | cut -f1))"
