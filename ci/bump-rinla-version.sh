#!/usr/bin/env bash
## Stamp rinla/DESCRIPTION's Version from the latest release tag.
##
## The R package version is maintained by hand, so it goes stale: devel
## carried 25.07.11.9000 for thirteen months while releases 26.08.07 and
## 26.08.22 shipped, which means anything installed from a checkout
## (R CMD INSTALL rinla, devtools::load_all) reported a version from the
## previous year.
##
## The release tags already are the source of truth and match the
## published packages one for one (tag v26.08.22 <-> INLA_26.08.22), so
## derive from them rather than keeping a second copy in step by hand:
##
##     Version: <latest tag without the leading v>.9000
##
## The .9000 suffix is R's convention for "development, after release X".
## It keeps a checkout build strictly newer than the release it follows,
## which is what makes an installed devel package win over the released
## one in version comparisons.
##
## Idempotent: writes only when the value actually changes, so running it
## on an already-current tree produces no diff. Use --check to report
## without writing (exit 1 if stale), which is what CI should call.
set -e -o pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DESC=$ROOT/rinla/DESCRIPTION
CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

[ -f "$DESC" ] || { echo "ERROR: no $DESC"; exit 1; }

## Nearest tag reachable from HEAD, not the lexically greatest: the tag
## list mixes zero-padded and unpadded forms (v26.08.20, v26.8.17) and a
## plain sort puts those in the wrong order.
TAG=$(git -C "$ROOT" describe --tags --abbrev=0 2>/dev/null || true)
[ -n "$TAG" ] || { echo "ERROR: no release tag reachable; cannot derive a version"; exit 1; }

## Strip whichever prefix the tag carries: Version_YY.MM.DD is the current
## form (matching the release/Version_* branches), v* the older one.
BASE=${TAG#Version_}
BASE=${BASE#v}
WANT="${BASE}.9000"
HAVE=$(awk -F': *' '/^Version:/ {print $2; exit}' "$DESC")

if [ "$HAVE" = "$WANT" ]; then
    echo "rinla/DESCRIPTION: Version $HAVE is current (tag $TAG)"
    exit 0
fi

if [ "$CHECK" = 1 ]; then
    echo "rinla/DESCRIPTION: Version $HAVE is stale, expected $WANT (tag $TAG)" >&2
    exit 1
fi

## In place, single line, and verified afterwards: a silent no-op here
## would leave the stale version in a package that looks freshly stamped.
tmp=$(mktemp)
awk -v want="$WANT" '/^Version:/ && !done { print "Version: " want; done=1; next } { print }' \
    "$DESC" > "$tmp"
mv "$tmp" "$DESC"

NOW=$(awk -F': *' '/^Version:/ {print $2; exit}' "$DESC")
[ "$NOW" = "$WANT" ] || { echo "ERROR: rewrite failed, Version is still $NOW"; exit 1; }
echo "rinla/DESCRIPTION: Version $HAVE -> $NOW (tag $TAG)"
