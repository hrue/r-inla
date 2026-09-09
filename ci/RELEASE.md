# Cutting a release

Everything below is enforced by the workflows. If a step is skipped the build
fails rather than publishing something wrong, so this is a description of what
the machinery already checks, not a list of manners.

## Before anything

**If libstiles changed, release sTiles FIRST.** Every lane fetches the solver
from a published sTiles release (`ci/fetch-stiles.sh`). Cutting INLA first
silently pairs the binaries with the previous solver, and a green board does
not reveal it. Check the "fetching ..." lines, or the `libstiles:` line in any
artifact's BUILDINFO.

## 1. Choose the version

    bash ci/bump-rinla-version.sh --release

Picks the next unused number for today and writes it to BOTH fields in
`rinla/DESCRIPTION`. It reads existing tags AND `rinla/NEWS.md`, because a
release branch can claim a number without ever tagging.

The version comes from a tag or from this command. Never from a clock: it does
not follow the commit date, the build date or the install date, so it cannot
move without a release and cannot go backwards.

`Version` and `Config/INLA/BinaryVersion` are always the same string. One
number identifies the R package and the binary that belongs with it.

**The rule this creates: every release must publish binaries.** The R package
asks for a binary release named by that field, so a tag without assets sends
users to a link that does not exist.

## 2. Write the NEWS entry

Add a `# INLA <version>` section at the top of `rinla/NEWS.md`. One short
bullet per user-visible change, in the style of the entries around it. This is
the only part nobody can generate: a changelog assembled from commit subjects
says what was edited, not what changed for a user.

A release with no entry fails the build.

## 3. Commit

    git add rinla/DESCRIPTION rinla/NEWS.md
    git commit -m "version <version>"
    git push origin feature/ci

`package/` regenerates itself in the commit hook. It is `rinla/` with the
symlinks resolved, and it is the only way a Windows user can install from the
repository.

## 4. Wait for green

Both workflows, on that commit. `build-inla` runs every platform and fits real
models on every shipped bundle: rgeneric and cgeneric, on both sparse-matrix
backends, on Windows, Linux x86_64, x86-64-v3, arm64 and macOS. `r-package`
builds the tarball and checks the metadata a release depends on.

Do not tag a red board. The release job will refuse anyway.

## 5. Tag

    git tag -a Version_<version> -m "Version_<version>"
    git push origin Version_<version>

The tag triggers a full rebuild on every platform and then publishes. It works
from any branch: the release job's only condition is that the ref is a tag.

## What the tag-time guards check

| guard | requires |
|---|---|
| Version is current | `DESCRIPTION` Version equals the tag |
| (same script) | `Config/INLA/BinaryVersion` equals Version |
| NEWS has an entry | a `# INLA <version>` section exists |
| package is in sync | `package/` matches `rinla/` |
| release needs | all eleven build and model-check jobs are green |

## After

Confirm the release actually carries assets, then confirm the pairing from the
artifacts rather than from intent: each one ships a BUILDINFO recording the
compiler, flags and the libstiles it was built against. That file and the CI
logs are the only records worth trusting.
