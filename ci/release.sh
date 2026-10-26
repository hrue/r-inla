#!/usr/bin/env bash
## Cut a release from the current branch in one command.
##
## Steps:
##   git checkout feature/ci && git pull --rebase origin feature/ci
##   bash ci/release.sh --dry-run          preview: version, NEWS, sTiles; changes nothing
##   bash ci/release.sh --news "Fixed X"   release; answer y once, then it runs alone:
##       commit version + NEWS, push, wait for r-package to go green, push the
##       tag, wait for build-inla (it publishes the release), check the assets,
##       copy the R package to r-inla/install. Ends with "done: <link>".
##   If libstiles changed, release sTiles first.
##   If devel has new commits: git merge origin/devel, then
##   bash ci/sync-package.sh && git add -A package && git commit -m "resync package"
##
## Options:
##   --dry-run, -n        show what would happen, change nothing
##   --news "text"        a NEWS line; repeat for more. Default: commit subjects
##                        since the last release
##   -y, --yes            no y/N prompt before the push
##   --no-mirror          leave r-inla/install alone
##   --mirror-only        only copy the newest release to r-inla/install
##   --no-wait-release    stop once the tag is pushed
##   -h, --help           this text
##
## If it stops:
##   r-package red     fix, commit, push, run again (same version)
##   build-inla red    gh run rerun <id> --failed (or fix and commit), run again
##   Ctrl-C, closed    run again; it resumes and never pushes anything twice
##
## On the release page: "Latest" is what inla.update() installs, a
## "Pre-release" reaches only inla.update(testing = TRUE).
##
## Needs gh (logged in, push rights) and rsync.
set -euo pipefail
## A CDPATH in the caller's shell makes cd print the directory, which corrupts
## every ROOT=$(cd ... && pwd) in the scripts this one calls.
unset CDPATH

YES=0
DRY=0
WAIT_RELEASE=1
REPO=""
NEWS=()
MIRROR=1
MIRROR_ONLY=0
MIRROR_REPO=${INLA_INSTALL_REPO:-r-inla/install}
## Files the install repository has and the package does not.
MIRROR_KEEP=(LICENSE README.md .gitignore)
## Assets the release job publishes: six bundles and the R tarball.
MIN_ASSETS=7

while [ $# -gt 0 ]; do
    case "$1" in
        -y|--yes)          YES=1 ;;
        -n|--dry-run)      DRY=1 ;;
        --no-wait-release) WAIT_RELEASE=0 ;;
        --repo)            shift; REPO=${1:-} ;;
        --news)            shift; NEWS+=("${1:-}") ;;
        --no-mirror)       MIRROR=0 ;;
        --mirror-only)     MIRROR_ONLY=1 ;;
        -h|--help)         awk 'NR > 1 && /^##/ { sub(/^## ?/, ""); print; next } NR > 1 { exit }' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

say()  { printf '==> %s\n' "$*"; }
## Release versions as pushed to origin. Local-only tags (a typo, a test) are
## not releases and must not count as the newest one.
release_versions() {
    git ls-remote --tags --refs origin 'Version_*' | sed 's#.*refs/tags/Version_##'
}
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
version() { sed -n 's/^Version:[[:space:]]*//p' rinla/DESCRIPTION | head -1; }
## Exact match: a section for 26.10.24-1 is not one for 26.10.24.
has_news() { grep -qE "^# INLA ${1//./\\.}([[:space:]]|\$)" rinla/NEWS.md; }
## Commit subjects since the last release, minus merges and housekeeping.
news_bullets() {
    if [ "${#NEWS[@]}" -gt 0 ]; then printf '%s\n' "${NEWS[@]}"; return; fi
    git log --no-merges --format='%s' ${NEWEST_TAG:+$NEWEST_TAG..}HEAD \
        | grep -viE '^(version |release|resync|rendered vignettes|bump|wip)|\[skip ci\]' \
        | awk 'NF && !seen[$0]++' | grep . || echo "Minor code improvements"
}

## Copy the package of release $1 into the install repository and push it.
## Built from the TAG, not the working tree, by the same ci/sync-package.sh
## that makes package/: that also brings the *.map demo data, which package/
## cannot track because the top-level .gitignore drops *.map. The rendered
## vignettes come from the "rendered vignettes for <tag>" commit the release
## job pushes after the tag; only rinla/inst/doc is taken from it, so nothing
## committed after the tag leaks in.
mirror_push() {
    local tag=$1 ver=${1#Version_} tmp vc url def
    tmp=$(mktemp -d); MTMP=$tmp
    trap 'git worktree remove --force "$MTMP/src" >/dev/null 2>&1 || true; rm -rf "$MTMP"' EXIT
    git worktree add -q --detach "$tmp/src" "$tag"
    git fetch -q origin
    vc=$(git log --all --format='%H %s' \
           | awk -v s="rendered vignettes for $tag [skip ci]" 'substr($0, 42) == s && !f { print $1; f = 1 }')
    if [ -n "$vc" ]; then
        git -C "$tmp/src" rm -rq --ignore-unmatch rinla/inst/doc
        git -C "$tmp/src" checkout -q "$vc" -- rinla/inst/doc
    else
        echo "WARNING: no rendered vignettes commit for $tag; the mirror gets the tag's inst/doc" >&2
    fi
    bash "$tmp/src/ci/sync-package.sh" >/dev/null
    [ "$(sed -n 's/^Version:[[:space:]]*//p' "$tmp/src/package/DESCRIPTION")" = "$ver" ] \
        || die "package/DESCRIPTION at $tag is not version $ver"

    ## Same protocol as origin, so the credentials that push the release
    ## also push here (an https URL prompts for a password over SSH setups).
    case "$(git remote get-url origin)" in
        git@*|ssh://*) url=git@github.com:$MIRROR_REPO.git ;;
        *)             url=https://github.com/$MIRROR_REPO.git ;;
    esac
    url=${INLA_INSTALL_URL:-$url}
    def=${INLA_INSTALL_BRANCH:-$(gh repo view "$MIRROR_REPO" --json defaultBranchRef --jq .defaultBranchRef.name)}
    git clone -q --depth 1 --branch "$def" "$url" "$tmp/mirror"
    local ex=(--exclude=/.git)
    for f in "${MIRROR_KEEP[@]}"; do ex+=("--exclude=/$f"); done
    rsync -a --delete "${ex[@]}" "$tmp/src/package/" "$tmp/mirror/"
    ## --force: the install repository's own .gitignore must not drop a file
    ## of the package.
    git -C "$tmp/mirror" add -A --force
    if git -C "$tmp/mirror" diff --cached --quiet; then
        say "$MIRROR_REPO already has $ver"
        return 0
    fi
    git -C "$tmp/mirror" commit -q -m "INLA $ver"
    git -C "$tmp/mirror" push -q origin "HEAD:$def"
    say "pushed INLA $ver to $MIRROR_REPO ($def)"
}

## Push access to the install repository; checked up front so a missing
## permission is a warning at the start, not a surprise at the end.
mirror_ok() {
    [ -n "${INLA_INSTALL_URL:-}" ] && return 0
    [ "$(gh api "repos/$MIRROR_REPO" --jq .permissions.push 2>/dev/null)" = true ]
}

cd "$(git rev-parse --show-toplevel)"
command -v gh >/dev/null || die "gh is not installed"
gh auth status >/dev/null 2>&1 || die "gh is not logged in (gh auth login)"

## owner/name from the origin URL, so this works from any clone or fork.
if [ -z "$REPO" ]; then
    REPO=$(git remote get-url origin | sed -E 's#^(git@[^:]+:|https?://[^/]+/)##; s#\.git$##')
fi
case "$REPO" in */*) ;; *) die "cannot tell the GitHub repository from origin; pass --repo owner/name" ;; esac

## Tracked changes would be swept into the release commit or block the
## checkout. Untracked files are fine.
git diff --quiet && git diff --cached --quiet \
    || die "the tree has uncommitted changes to tracked files; commit or stash them"

if [ "$MIRROR_ONLY" = 1 ]; then
    git fetch -q --tags origin
    T=$(release_versions | sort -V | tail -1)
    [ -n "$T" ] || die "no release tag"
    gh release view "Version_$T" -R "$REPO" >/dev/null 2>&1 || die "Version_$T has no published release"
    mirror_ok || die "no push access to $MIRROR_REPO"
    [ "$DRY" = 1 ] && { say "would copy Version_$T to $MIRROR_REPO"; exit 0; }
    mirror_push "Version_$T"
    exit 0
fi

## Release the branch we are on. It has to be one r-package runs on, or there
## is no green run to wait for.
BRANCH=$(git symbolic-ref -q --short HEAD) || die "detached HEAD; check out the branch to release"
RUNS_ON=$(sed -n 's/^ *branches: *\[\(.*\)\].*/\1/p' .github/workflows/r-package.yml | head -1 | tr -d ' ')
case ",$RUNS_ON," in
    *",$BRANCH,"*) ;;
    *) die "r-package does not run on $BRANCH (only $RUNS_ON), so nothing would go green" ;;
esac

say "repository $REPO, branch $BRANCH"
git fetch -q --tags origin
git merge -q --ff-only "origin/$BRANCH" || die "$BRANCH has diverged from origin/$BRANCH"

if [ "$MIRROR" = 1 ] && ! mirror_ok; then
    echo "WARNING: no push access to $MIRROR_REPO; the release will not reach it" >&2
    MIRROR=0
fi

## Tagged at the tip of the branch with no release yet: the tag build failed
## or is still running. Resume there instead of cutting another version.
V=$(version)
NEWEST=$(release_versions | sort -V | tail -1)
NEWEST_TAG=${NEWEST:+Version_$NEWEST}
RESUME=0
if [ "$(git rev-parse -q --verify "refs/tags/Version_$V^{commit}" || true)" = "$(git rev-parse HEAD)" ] \
   && ! gh release view "Version_$V" -R "$REPO" >/dev/null 2>&1; then
    RESUME=1
    say "Version_$V is tagged here but not published; resuming at the release build"
fi

## Releasing a branch that lacks the newest release would publish older code
## as "latest". And with nothing new since it there is nothing to release.
if [ "$RESUME" = 0 ] && [ -n "$NEWEST_TAG" ]; then
    git merge-base --is-ancestor "$NEWEST_TAG" HEAD \
        || die "$BRANCH does not contain $NEWEST_TAG; merge the branch it was cut from first"
    [ "$(git rev-list --count "$NEWEST_TAG..HEAD")" -gt 0 ] \
        || die "nothing on $BRANCH since $NEWEST_TAG"
fi

## A version is spent if it is tagged, or if it is not newer than the newest
## release (a tag cut from another branch, or a DESCRIPTION left behind its
## tag). Otherwise it is this release, typically from an interrupted run.
if [ "$RESUME" = 1 ]; then
    :
elif git rev-parse -q --verify "refs/tags/Version_$V" >/dev/null \
   || [ "$(printf '%s\n%s\n' "$V" "$NEWEST" | sort -V | tail -1)" != "$V" ]; then
    if [ "$DRY" = 1 ]; then
        say "$V is not newer than Version_$NEWEST; would run ci/bump-rinla-version.sh --release"
        V="<next free for $(date +%y.%m.%d)>"
    else
        bash ci/bump-rinla-version.sh --release
        V=$(version)
    fi
else
    say "keeping $V, it is not released yet"
fi
TAG="Version_$V"

## The sTiles solver a tag build will pair with. Release sTiles first if
## libstiles changed, or the bundles get the previous one.
## Same repository the build fetches from (ci/fetch-stiles.sh), so the two
## cannot disagree.
STILES_REPO=${STILES_REPO:-$(sed -n 's/^STILES_REPO=\${STILES_REPO:-\(.*\)}$/\1/p' ci/fetch-stiles.sh)}
ST=$(gh api "repos/$STILES_REPO/releases/latest" --jq '.tag_name + " (" + (.published_at[:10]) + ")"' 2>/dev/null || echo unknown)
say "bundles will be built against sTiles $ST"

if [ "$DRY" = 1 ] && [ "$RESUME" = 1 ]; then
    say "would wait for the build-inla run of $TAG and check the release"
    exit 0
elif [ "$DRY" = 1 ]; then
    has_news "$V" \
        && say "NEWS already has a section for $V" \
        || { say "would add to rinla/NEWS.md:"; news_bullets | sed 's/^/      * /'; }
    say "would commit, push $BRANCH, wait for r-package, then push tag $TAG"
    [ "$MIRROR" = 1 ] && say "then copy the released package to $MIRROR_REPO"
    last=$(gh run list -R "$REPO" --workflow r-package.yml --branch "$BRANCH" --limit 1 \
             --json conclusion,headSha --jq '.[0] | (.headSha[:9]) + " " + .conclusion' 2>/dev/null || true)
    say "last r-package run on $BRANCH: ${last:-none}"
    exit 0
fi

## Wait for the run of THIS commit on THIS ref. "The latest run" can still be
## the previous commit's for a few seconds after a push, and one commit has a
## run per ref: the branch push and the tag push are separate runs.
wait_run() {   # wait_run <workflow> <sha> <branch or tag> <minutes to appear> <what to do if red>
    local wf=$1 sha=$2 ref=$3 tries=$(( $4 * 4 )) hint=$5 id=""
    while [ "$tries" -gt 0 ]; do
        id=$(gh run list -R "$REPO" --workflow "$wf" --commit "$sha" --branch "$ref" --limit 1 \
               --json databaseId --jq '.[0].databaseId // empty' 2>/dev/null || true)
        [ -n "$id" ] && break
        tries=$((tries - 1))
        sleep 15
    done
    [ -n "$id" ] || die "no $wf run appeared for $ref at ${sha:0:9}"
    say "watching $wf run $id (https://github.com/$REPO/actions/runs/$id)"
    gh run watch "$id" -R "$REPO" --exit-status --interval 30 >/dev/null \
        || die "$wf is red for ${sha:0:9}: $hint"
    say "$wf is green"
}

SHA=$(git rev-parse HEAD)
if [ "$RESUME" = 0 ]; then

if ! has_news "$V"; then
    { printf '# INLA %s\n' "$V"; news_bullets | sed 's/^/* /'; echo; cat rinla/NEWS.md; } > rinla/NEWS.md.tmp
    mv rinla/NEWS.md.tmp rinla/NEWS.md
    say "NEWS for $V:"; sed -n '2,/^$/p' rinla/NEWS.md | sed '/^$/d; s/^/      /'
fi

## package/ in the same commit. The commit hook normally does this, but not
## when core.hooksPath points elsewhere, and the r-package guard fails on a
## stale package/. The sync reads tracked files, so stage rinla/ first.
git add rinla/DESCRIPTION rinla/NEWS.md
bash ci/sync-package.sh >/dev/null
git add -A package
git diff --cached --quiet || git commit -q -m "Version $V"

SHA=$(git rev-parse HEAD)
AHEAD=$(git rev-list --count "origin/$BRANCH..HEAD")
if [ "$AHEAD" -gt 0 ]; then
    echo
    git --no-pager log --oneline "origin/$BRANCH..HEAD"
    echo
    if [ "$YES" != 1 ]; then
        [ -t 0 ] || die "not a terminal: re-run with -y to push without asking"
        read -r -p "Push these to $REPO $BRANCH and release $TAG? [y/N] " ans
        case "$ans" in y|Y|yes) ;; *) die "stopped; the commits are local only" ;; esac
    fi
    git push -q origin "$BRANCH"
    say "pushed $BRANCH at ${SHA:0:9}"
else
    say "$BRANCH is already pushed at ${SHA:0:9}"
fi

wait_run r-package.yml "$SHA" "$BRANCH" 10 "fix it, push, and run this again"

## Tag. An existing tag is fine only if it is this commit.
if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
    [ "$(git rev-parse "$TAG^{commit}")" = "$SHA" ] || die "$TAG already exists on another commit"
else
    git tag -a "$TAG" -m "$TAG"
fi
if git ls-remote --exit-code --tags origin "refs/tags/$TAG" >/dev/null 2>&1; then
    say "$TAG is already on origin"
else
    git push -q origin "$TAG"
    say "pushed $TAG"
fi

fi   # RESUME = 0

if [ "$WAIT_RELEASE" = 0 ]; then
    say "done; the release job publishes when build-inla finishes"
    [ "$MIRROR" = 1 ] && say "after that: bash ci/release.sh --mirror-only"
    exit 0
fi

## The tag starts the full board, which ends in the release job. The tag
## points at the commit, so the run is found by the same sha.
wait_run build-inla.yml "$SHA" "$TAG" 10 \
    "re-run the failed jobs (gh run rerun <id> --failed -R $REPO), then run this again"

n=$(gh release view "$TAG" -R "$REPO" --json assets --jq '.assets | length' 2>/dev/null || echo 0)
[ "$n" -ge "$MIN_ASSETS" ] || die "release $TAG has $n assets, expected $MIN_ASSETS"
rows=$(gh release view "$TAG" -R "$REPO" --json body --jq .body | grep -c '^| inla' || true)
latest=$(gh api "repos/$REPO/releases/latest" --jq .tag_name 2>/dev/null || true)
say "release $TAG: $n assets, $rows rows in the notes table, latest = $latest"
[ "$rows" -gt 0 ] || echo "WARNING: the notes table is empty; check the manifest artifacts" >&2
[ "$latest" = "$TAG" ] || echo "WARNING: releases/latest is $latest, not $TAG" >&2
if [ "$MIRROR" = 1 ]; then
    mirror_push "$TAG"
fi
say "done: https://github.com/$REPO/releases/tag/$TAG"
