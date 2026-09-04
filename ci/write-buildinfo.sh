#!/usr/bin/env bash
## Write a BUILDINFO file that travels inside the artifact: the exact
## compiler, flags, libraries and toolchain settings this binary was built
## with. The point is that the answer to "how was this compiled?" ships
## with the binary instead of living in CI logs that expire.
##
##   ci/write-buildinfo.sh <outfile> <CC> <FLAGS> <BLAS-description>
set -e

OUT=$1; CC=$2; FLAGS=$3; BLAS=$4
ROOT=$(cd "$(dirname "$0")/.." && pwd)

{
    echo "== inla build =="
    echo "date:       $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "commit:     $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "compiler:   $($CC --version 2>/dev/null | head -1)"
    echo "flags:      $FLAGS"
    echo "blas:       $BLAS"
    echo "libstiles:  $(cat "$ROOT"/local*/stiles/RELEASE 2>/dev/null | head -1 || echo none)"
    echo "os:         $(uname -sm) $( (. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME") || sw_vers -productVersion 2>/dev/null || true)"
    echo
    echo "== shared toolchain file (ci/toolchain.env) =="
    grep -vE '^\s*#|^\s*$' "$ROOT/ci/toolchain.env" 2>/dev/null || echo "(absent)"
} > "$OUT"
echo "buildinfo: $OUT"
