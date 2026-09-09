#!/usr/bin/env bash
## Keep rinla/DESCRIPTION's Version, and the binary release it names, correct.
##
## A version comes from exactly two places: a release tag, or --release when a
## release is being cut. Never from a clock. It used to be stamped with the
## HEAD commit's date, or today's date from the commit hook, which meant it
## moved on every commit with no release happening, two different builds could
## claim the same released number, and it could go BACKWARDS when a machine's
## date disagreed, so a correct dependency like INLA (>= 26.09.08) was
## rejected for a package that really was 26.09.08.
##
##     at a tag        the tag's version
##     --release       the next unused number for today
##     anything else   unchanged: it belongs to the last release
##
## Version and Config/INLA/BinaryVersion are always the same string, so one
## number identifies the R package and the binary that belongs with it. The
## rule that creates: every release must publish binaries, because that field
## names a release tag which has to exist.
##
## No development suffix: the built binary reports this exact string, so
## anything like .9000 would show up in `inla -V`. The consequence, worth
## knowing: a package installed from a checkout carries the LAST RELEASE's
## version, so R cannot tell it apart from the released one by version alone.
## Use the commit, not the version, to identify a checkout build.
##
## Idempotent: writes only when the value actually changes, so running it on an
## already-current tree produces no diff. --check reports without writing and
## exits 1 if wrong, which is what CI calls. --check and --release are mutually
## exclusive: one verifies the number in use, the other picks an unused one.
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
        ## Pick the next FREE version for a release today: the plain date if it
        ## has never been used, else <date>-N with the lowest N that is free.
        ## Without this the suffix was hand-edited, which is how a number that a
        ## release branch had already claimed could be reused.
        --release) RELEASE=1 ;;
        ## Also stamp Config/INLA/BinaryVersion, the field that says which
        ## BINARY release this R package needs. It is deliberately NOT tied to
        ## Version: R-only edits are frequent and need no new solver. Pass this
        ## exactly when the C sources changed, which is when a new binary is
        ## genuinely required; the pre-commit hook decides that by looking at
        ## what is staged.
        ## Accepted and ignored: BinaryVersion now always tracks Version,
        ## so there is nothing to opt into. Kept so an existing caller
        ## does not fail on an unknown option.
        --binary) : ;;
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
## --release and --check contradict each other: --release picks the next
## UNUSED number, so it always reports the current one as stale. The guard is
## plain --check, which verifies the version that is already set.
if [ "$RELEASE" = 1 ] && [ "$CHECK" = 1 ]; then
    echo "usage: --release and --check are mutually exclusive" >&2
    exit 2
fi

## Read the version that is already set BEFORE deciding what it should be:
## when nothing decides otherwise, the answer is that it does not change.
HAVE=$(awk -F': *' '/^Version:/ {print $2; exit}' "$DESC")

if [ "$RELEASE" = 1 ]; then
    ## Both TAGS and NEWS.md are consulted. Tags alone are not enough: a
    ## release branch can claim a number and be merged without ever tagging,
    ## which is exactly how 26.09.08-1 and -2 came to exist with no tag.
    _base=$(date +%y.%m.%d)
    _used=$( { git -C "$ROOT" tag --list "Version_${_base}*" | sed 's/^Version_//'
               grep -oE "^# INLA ${_base}(-[0-9]+)?" "$ROOT/rinla/NEWS.md" 2>/dev/null \
                   | sed 's/^# INLA //'; } | sort -u )
    if ! printf '%s\n' "$_used" | grep -qx "$_base"; then
        WANT="$_base"
    else
        _n=1
        while printf '%s\n' "$_used" | grep -qx "$_base-$_n"; do
            _n=$((_n + 1))
        done
        WANT="$_base-$_n"
    fi
elif [ -n "$(git -C "$ROOT" tag --points-at HEAD 2>/dev/null)" ]; then
    ## At a tag the version IS that release.
    WANT="$BASE"
else
    ## Not at a tag: the version does NOT change. It belongs to the last
    ## release and stays there until the next one is cut.
    ##
    ## It used to come from the clock here (today's date via --now from the
    ## commit hook, or the HEAD commit's date otherwise), so it moved on every
    ## commit with no release happening. Two different versions could then
    ## describe the same released code, and the number could go BACKWARDS when
    ## a machine's date disagreed, so Depends: INLA (>= 26.09.08) could be
    ## rejected for a package that really was 26.09.08. A version now comes
    ## from a tag or from --release, never from a clock.
    WANT="$HAVE"
fi
[ -n "$WANT" ] || { echo "ERROR: could not derive a version"; exit 1; }

## Config/INLA/BinaryVersion is now ALWAYS the same string as Version.
##
## It used to move only when the C sources changed, so that an R-only fix did
## not force a binary release. The cost was two similar-looking dates that
## disagreed (Version 26.09.07-1 against BinaryVersion 26.09.07) and nobody
## could say which number identified what they had. One number now identifies
## the pair: the R package and the binary it belongs with.
##
## THE RULE THIS CREATES: every release must publish binaries. The R package
## asks for a binary release named by this field, so a release that bumps it
## without publishing the assets sends users to a tag that does not exist.
##
## Takes the version as an argument rather than reading $WANT, because a
## same-day re-release keeps its suffix (26.09.07-1) and the binary field has
## to carry the suffix too, not the bare date.
stamp_binary() {
    _want=$1
    bhave=$(awk -F': *' '/^Config\/INLA\/BinaryVersion:/ {print $2; exit}' "$DESC")
    [ -n "$bhave" ] || return 0
    [ "$bhave" = "$_want" ] && return 0
    if [ "$CHECK" = 1 ]; then
        echo "rinla/DESCRIPTION: Config/INLA/BinaryVersion is $bhave, expected $_want" >&2
        exit 1
    fi
    btmp=$(mktemp)
    awk -v want="$_want" '/^Config\/INLA\/BinaryVersion:/ && !d { print "Config/INLA/BinaryVersion: " want; d=1; next } { print }' \
        "$DESC" > "$btmp"
    mv "$btmp" "$DESC"
    echo "rinla/DESCRIPTION: Config/INLA/BinaryVersion $bhave -> $_want"
}

## A SECOND release on the same day is written <date>-N (26.09.06-2), to match
## a Version_26.09.06-2 tag: the date alone is taken. Such a version is
## deliberate and current, so leave it alone. Without this the hook rewrote it
## back to the plain date on the next commit, and the guard called it stale on
## every branch push, so a same-day re-release could not be held.
case "$HAVE" in
    "$WANT"-[0-9]*)
        echo "rinla/DESCRIPTION: Version $HAVE is current (same-day re-release of $WANT)"
        ## the suffixed version is the effective one, so the binary field takes it
        stamp_binary "$HAVE"
        exit 0
        ;;
esac

if [ "$HAVE" = "$WANT" ]; then
    echo "rinla/DESCRIPTION: Version $HAVE is current (tag $TAG)"
    stamp_binary "$HAVE"
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
stamp_binary "$WANT"
