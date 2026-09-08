#!/usr/bin/env bash
## Regenerate package/ from rinla/.
##
## package/ is rinla/ with every symlink resolved into the file it points to.
## It exists because Windows cannot extract this repository: tar.exe fails on
## each of the 211 symlinks with "Can't create ...: Invalid argument" and then
## aborts, so remotes::install_github() dies before it reaches any subdir.
## .gitattributes keeps the symlink-bearing paths out of the generated
## archives, and this directory carries the package instead:
##
##     remotes::install_github("hrue/r-inla", subdir = "package")
##
## Never edit package/ by hand. Edit rinla/, run this, commit both. The
## r-package workflow runs it with --check and fails if they have drifted,
## because a copy that silently falls behind is worse than no copy: it would
## ship an old package to exactly the users who cannot install any other way.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)

## The file list comes from GIT, not from the working tree.
##
## Copying rinla/ wholesale mirrors whatever happens to be on that machine,
## which includes build artifacts nobody tracks: NAMESPACE.backup, R/TAGS,
## regenerated vignette figures. Those got copied into package/ and committed
## from one developer's tree, and every clean checkout then reported package/
## as out of date by exactly those files. The shipped package must not depend
## on whose laptop ran the sync.
##
## Rebuilt from empty each time rather than synced in place: it makes the
## result a pure function of the tracked tree, and removes any question about
## rsync's comparison rules (the previous version needed --checksum because a
## fresh checkout gives every file the same mtime).
rm -rf "$ROOT/package"
mkdir -p "$ROOT/package"
git -C "$ROOT/rinla" ls-files -z \
    | rsync -a --copy-links --files-from=- --from0 "$ROOT/rinla/" "$ROOT/package/"

## No symlink may survive: --copy-links resolves them, but a link whose target
## is missing is skipped rather than resolved, which would put the problem
## straight back into the archive.
if find "$ROOT/package" -type l | grep -q .; then
    echo "ERROR: package/ still contains symlinks (broken targets in rinla/?):" >&2
    find "$ROOT/package" -type l | sed 's/^/    /' >&2
    exit 1
fi

if [ "${1:-}" = "--check" ]; then
    if [ -n "$(git -C "$ROOT" status --porcelain -- package)" ]; then
        echo "::error::package/ is out of date. Run 'bash ci/sync-package.sh' and commit the result."
        git -C "$ROOT" status --porcelain -- package | head -20
        exit 1
    fi
    echo "package/ is in sync with rinla/"
else
    echo "package/ regenerated from rinla/"
fi
