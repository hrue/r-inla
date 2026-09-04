#!/usr/bin/env bash
## Draft a NEWS.md section for the version being released.
##
## NEWS entries are summaries a reader can act on ("Minor code improvements",
## "Last build on R-4.6"), not a commit log. So this DRAFTS, it does not
## decide: it collects the commit subjects since the last release, drops the
## ones that are never news (merges, duplicates, CI-only churn), and leaves a
## section to edit down. The judgement of what mattered stays with the person
## writing it; what is automated is the part worth automating, which is not
## starting from an empty file and not missing the section entirely.
##
##   ci/news-draft.sh              print the draft
##   ci/news-draft.sh --insert     prepend it to rinla/NEWS.md, then edit
##   ci/news-draft.sh --since v26.08.20   from a specific tag
##
## The release guard in the r-package workflow requires a "# INLA <version>"
## section for the version in DESCRIPTION, so run this before tagging.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
NEWS=$ROOT/rinla/NEWS.md
DESC=$ROOT/rinla/DESCRIPTION

INSERT=0
SINCE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --insert) INSERT=1 ;;
        --since)  shift; SINCE=${1:-} ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

VER=$(sed -n 's/^Version:[[:space:]]*//p' "$DESC" | head -1)
[ -n "$VER" ] || { echo "ERROR: no Version in $DESC" >&2; exit 1; }

if [ -z "$SINCE" ]; then
    SINCE=$(git -C "$ROOT" describe --tags --abbrev=0 2>/dev/null || true)
fi
RANGE=${SINCE:+$SINCE..}HEAD

## Already there? Then say so rather than stacking a second section for the
## same version, which is the obvious way to make NEWS.md wrong.
if grep -qE "^#[[:space:]]*INLA[[:space:]]+${VER}([[:space:]]|\$)" "$NEWS"; then
    echo "rinla/NEWS.md already has a section for $VER; nothing to draft."
    exit 0
fi

## What is never news:
##  - merge commits (no content of their own)
##  - the release commits themselves
##  - anything marked as skipping CI, which is by definition plumbing
## Duplicates collapse: the same subject often appears from a merge of the
## same work through two branches.
subjects=$(git -C "$ROOT" log --no-merges --format='%s' $RANGE 2>/dev/null \
    | grep -vEi '^(Release|Merge |bump|wip\b)' \
    | sed 's/[[:space:]]*\[skip ci\][[:space:]]*//I' \
    | grep -vE '^[[:space:]]*$' \
    | awk '!seen[$0]++')

n=$(printf '%s\n' "$subjects" | grep -c . || true)

{
    echo "# INLA $VER"
    if [ "$n" -eq 0 ]; then
        echo "* "
    else
        printf '%s\n' "$subjects" | sed 's/^/* /'
    fi
} > /tmp/news-draft.$$

if [ "$INSERT" = 1 ]; then
    ## Prepend, keeping a blank line between sections, and never in place:
    ## write beside and rename, so an interrupted run cannot truncate NEWS.md.
    tmp=$(mktemp)
    { cat /tmp/news-draft.$$; echo; cat "$NEWS"; } > "$tmp"
    mv "$tmp" "$NEWS"
    rm -f /tmp/news-draft.$$
    echo "Prepended a draft section for $VER to rinla/NEWS.md ($n entries from ${SINCE:-the start})."
    echo "EDIT IT before committing: these are commit subjects, not release notes."
else
    cat /tmp/news-draft.$$
    rm -f /tmp/news-draft.$$
    echo
    echo "## $n commits since ${SINCE:-the start}. Re-run with --insert to prepend this to rinla/NEWS.md."
fi
