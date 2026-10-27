## Sourced by the ci/deps-*.sh scripts and ci/fetch-stiles.sh.
##
## Every lane downloads hundreds of packages and a dozen source trees, and
## each of those servers has bad minutes: a crawling Ubuntu mirror, a 504 from
## a GitHub release download, GitLab "unable to handle this request due to
## load". One failed download used to fail a twenty minute job. After this,
## the downloads in the sourcing script retry instead.

## Root in the containers, sudo on the hosted runners.
_SUDO=
[ "$(id -u)" = 0 ] || _SUDO=sudo

## curl and wget read these files themselves, so every call retries without
## being edited. 504, 429 and dropped connections are all retried. A private
## directory, so running a deps script by hand never touches ~/.curlrc.
_NETRC_DIR=$(mktemp -d)
export CURL_HOME=$_NETRC_DIR WGETRC=$_NETRC_DIR/.wgetrc
cat > "$CURL_HOME/.curlrc" <<'EOF'
retry = 6
retry-delay = 15
retry-max-time = 600
connect-timeout = 30
EOF
cat > "$WGETRC" <<'EOF'
tries = 6
waitretry = 15
retry_connrefused = on
timeout = 60
EOF

## System package managers only on a CI runner, never on someone's machine.
if [ "${GITHUB_ACTIONS:-}" = true ]; then
    ## apt: retry each file, and give up on a connection that stalls.
    if [ -d /etc/apt/apt.conf.d ]; then
        printf 'Acquire::Retries "5";\nAcquire::http::Timeout "60";\nAcquire::https::Timeout "60";\n' \
            | $_SUDO tee /etc/apt/apt.conf.d/80-ci-retries >/dev/null
    fi
    ## dnf: retry, and leave a mirror that crawls below 50 kB/s for 60 s.
    for _conf in /etc/dnf/dnf.conf /etc/yum.conf; do
        if [ -f "$_conf" ] && ! grep -q '^minrate=' "$_conf"; then
            printf 'retries=10\ntimeout=60\nminrate=50k\n' | $_SUDO tee -a "$_conf" >/dev/null
        fi
    done
fi

## apt-get with a time limit per attempt. Ubuntu's Azure mirror sometimes
## serves at a crawl without ever failing a file, which apt cannot detect;
## an install that normally takes two minutes then ate the whole step. So
## stop after seven minutes and go again, from the main archive.
apt_get() {
    local i
    for i in 1 2 3; do
        $_SUDO timeout 420 apt-get "$@" && return 0
        echo "apt-get $1 failed or stalled (attempt $i/3)" >&2
        [ "$i" = 3 ] && return 1
        $_SUDO dpkg --configure -a || true
        [ "${GITHUB_ACTIONS:-}" = true ] && $_SUDO sed -i \
            's#azure\.archive\.ubuntu\.com#archive.ubuntu.com#g' \
            /etc/apt/sources.list /etc/apt/sources.list.d/*.list \
            /etc/apt/sources.list.d/*.sources 2>/dev/null || true
        $_SUDO timeout 300 apt-get update || true
    done
}

## git clone with retries, clearing the half made target in between.
## Shadows git only in the shell that sources this file.
git() {
    if [ "${1:-}" != clone ]; then
        command git "$@"
        return
    fi
    local i dest=${!#}
    for i in 1 2 3 4 5; do
        command git "$@" && return 0
        [ "$i" = 5 ] && return 1
        echo "git clone failed (attempt $i/5), again in $((i * 20)) s" >&2
        case "$dest" in *://*|git@*) ;; *) rm -rf "$dest" ;; esac
        sleep $((i * 20))
    done
}

## GKlib for METIS, fetched with git rather than as the anonymous tarball
## METIS's CMake downloads. Anonymous downloads are throttled per runner IP,
## while the git fetch goes through the job token (see the insteadOf in the
## workflow). Prints the CMake flag that hands it over; empty if not needed.
gklib_for_metis() {
    local metis=$1 sha
    sha=$(grep -oE 'GKlib/archive/[0-9a-f]{40}' "$metis/CMakeLists.txt" 2>/dev/null \
          | head -1 | cut -d/ -f3)
    [ -n "$sha" ] || return 0
    rm -rf /tmp/gklib
    command git init -q /tmp/gklib
    local i
    for i in 1 2 3 4 5; do
        command git -C /tmp/gklib fetch -q --depth 1 https://github.com/scivision/GKlib "$sha" \
            && command git -C /tmp/gklib checkout -q FETCH_HEAD \
            && { echo "-DFETCHCONTENT_SOURCE_DIR_GKLIB=/tmp/gklib"; return 0; }
        echo "GKlib fetch failed (attempt $i/5)" >&2
        sleep $((i * 20))
    done
    return 0   # fall back to METIS's own download
}
