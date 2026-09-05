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
##     Version: <YY.MM.DD of HEAD, or the tag when HEAD is tagged>
##
## No development suffix: the built binary reports this exact string.
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
NOW=0
BINARY=0
for a in "$@"; do
    case "$a" in
        --check) CHECK=1 ;;
        ## Date the version from TODAY rather than from HEAD. For the
        ## pre-commit hook: at that moment HEAD is still the PREVIOUS commit,
        ## so the HEAD-based date is one commit behind and the first commit of
        ## a new day would ship a stale version, which is exactly the case the
        ## CI guard keeps catching.
        --now)   NOW=1 ;;
        ## Also stamp Config/INLA/BinaryVersion, the field that says which
        ## BINARY release this R package needs. It is deliberately NOT tied to
        ## Version: R-only edits are frequent and need no new solver. Pass this
        ## exactly when the C sources changed, which is when a new binary is
        ## genuinely required; the pre-commit hook decides that by looking at
        ## what is staged.
        --binary) BINARY=1 ;;
        *) echo "usage: $0 [--check] [--now]" >&2; exit 2 ;;
    esac
done

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
## Exactly at a tag -> that release's version. Otherwise the version is the
## HEAD commit's DATE, in the same YY.MM.DD form the releases use. There is
## deliberately NO ".9000" development suffix: the build scripts read this
## same string out of DESCRIPTION and compile it into the binary, so it is
## what a user sees from BOTH packageVersion("INLA") and `inla -V`, and a
## suffix there reads as noise. The date always sorts above the last tag, so
## R still sees an upgrade, and it never claims to BE a release the way a
## bare tag on a later commit would.
if [ "$NOW" = 1 ]; then
    WANT=$(date +%y.%m.%d)
elif [ -n "$(git -C "$ROOT" tag --points-at HEAD 2>/dev/null)" ]; then
    WANT="$BASE"
else
    WANT=$(git -C "$ROOT" log -1 --format=%cd --date=format:%y.%m.%d 2>/dev/null)
fi
[ -n "$WANT" ] || { echo "ERROR: could not derive a version"; exit 1; }
HAVE=$(awk -F': *' '/^Version:/ {print $2; exit}' "$DESC")


## Config/INLA/BinaryVersion: which BINARY release this R package needs. Only
## touched with --binary, i.e. when the C sources changed in this commit.
stamp_binary() {
    [ "$BINARY" = 1 ] || return 0
    bhave=$(awk -F': *' '/^Config\/INLA\/BinaryVersion:/ {print $2; exit}' "$DESC")
    [ -n "$bhave" ] || return 0
    [ "$bhave" = "$WANT" ] && return 0
    btmp=$(mktemp)
    awk -v want="$WANT" '/^Config\/INLA\/BinaryVersion:/ && !d { print "Config/INLA/BinaryVersion: " want; d=1; next } { print }' \
        "$DESC" > "$btmp"
    mv "$btmp" "$DESC"
    echo "rinla/DESCRIPTION: Config/INLA/BinaryVersion $bhave -> $WANT (C sources changed)"
}

if [ "$HAVE" = "$WANT" ]; then
    echo "rinla/DESCRIPTION: Version $HAVE is current (tag $TAG)"
    stamp_binary
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

NOW_HAVE=$(awk -F': *' '/^Version:/ {print $2; exit}' "$DESC")
[ "$NOW_HAVE" = "$WANT" ] || { echo "ERROR: rewrite failed, Version is still $NOW_HAVE"; exit 1; }
echo "rinla/DESCRIPTION: Version $HAVE -> $WANT (tag $TAG)"
stamp_binary
