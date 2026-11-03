#!/usr/bin/env bash
## Full build of a branch on GitHub, in one command. Publishes nothing.
##
##   bash ci/check-build.sh                 build the branch you are on, wait
##   bash ci/check-build.sh devel           build devel
##   bash ci/check-build.sh --vignettes     also render the vignettes (r-package)
##   bash ci/check-build.sh --no-wait       start it and print the link
##
## It runs build-inla by hand on the branch as it is on GitHub: all six
## bundles, a model on Linux, macOS and Windows, and the distro tests. No tag,
## no release, no change to any branch. Ends with "green" or the failed jobs.
## About 20 minutes; --vignettes adds r-package, about 15 more, in parallel.
##
## Needs gh (logged in).
set -euo pipefail

REPO=""
WAIT=1
VIGNETTES=0
BRANCH=""
while [ $# -gt 0 ]; do
    case "$1" in
        --no-wait)    WAIT=0 ;;
        --vignettes)  VIGNETTES=1 ;;
        --repo)       shift; REPO=${1:-} ;;
        -h|--help)    sed -n '2,14p' "$0"; exit 0 ;;
        -*)           echo "unknown option: $1" >&2; exit 2 ;;
        *)            BRANCH=$1 ;;
    esac
    shift
done

say() { printf '==> %s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

cd "$(git rev-parse --show-toplevel)"
command -v gh >/dev/null || die "gh is not installed"
gh auth status >/dev/null 2>&1 || die "gh is not logged in (gh auth login)"
if [ -z "$REPO" ]; then
    REPO=$(git remote get-url origin | sed -E 's#^(git@[^:]+:|https?://[^/]+/)##; s#\.git$##')
fi
case "$REPO" in */*) ;; *) die "cannot tell the GitHub repository from origin; pass --repo owner/name" ;; esac

[ -n "$BRANCH" ] || BRANCH=$(git symbolic-ref -q --short HEAD) || die "detached HEAD; name the branch"

## GitHub builds what is pushed, not this tree. Say so when they differ.
SHA=$(gh api "repos/$REPO/commits/$BRANCH" --jq .sha 2>/dev/null) || die "no branch $BRANCH on $REPO"
say "building $BRANCH at ${SHA:0:9} on $REPO"
if [ "$(git symbolic-ref -q --short HEAD || true)" = "$BRANCH" ] && [ "$(git rev-parse HEAD)" != "$SHA" ]; then
    say "NOTE: your local $BRANCH is at $(git rev-parse --short HEAD); push first to build that"
fi

## Newest manual run of a workflow for this commit, so an older run on the
## same commit is not mistaken for the one just started.
newest() {
    gh run list -R "$REPO" --workflow "$1" --event workflow_dispatch --branch "$BRANCH" \
        --commit "$SHA" --limit 1 --json databaseId --jq '.[0].databaseId // 0' 2>/dev/null || echo 0
}
start() {   # start <workflow> [gh workflow run args]
    local wf=$1 before id k; shift
    before=$(newest "$wf")
    gh workflow run "$wf" -R "$REPO" --ref "$BRANCH" "$@" >/dev/null
    for k in $(seq 1 40); do
        sleep 3
        id=$(newest "$wf")
        [ "$id" -gt "$before" ] 2>/dev/null && { echo "$id"; return 0; }
    done
    die "$wf was started but no run appeared"
}

B=$(start build-inla.yml)
say "build-inla: https://github.com/$REPO/actions/runs/$B"
R=""
if [ "$VIGNETTES" = 1 ]; then
    R=$(start r-package.yml -f vignettes=true)
    say "r-package (vignettes): https://github.com/$REPO/actions/runs/$R"
fi
[ "$WAIT" = 1 ] || exit 0

## Red: list the failed jobs with their links, so the next step is one click.
report() {   # report <run id> <name>
    if gh run watch "$1" -R "$REPO" --exit-status --interval 30 >/dev/null 2>&1; then
        say "$2 is green"
        return 0
    fi
    say "$2 is RED: https://github.com/$REPO/actions/runs/$1"
    gh run view "$1" -R "$REPO" --json jobs \
        --jq '.jobs[] | select(.conclusion == "failure") | "      \(.name)\n        \(.url)"'
    return 1
}

say "waiting (about 20 minutes; Ctrl-C stops waiting, not the build)"
ok=0
[ -z "$R" ] || report "$R" "r-package (vignettes)" || ok=1
report "$B" build-inla || ok=1
[ "$ok" = 0 ] && say "all green for $BRANCH at ${SHA:0:9}"
exit "$ok"
