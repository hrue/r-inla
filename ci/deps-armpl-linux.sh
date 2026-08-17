#!/usr/bin/env bash
## Install Arm Performance Libraries (aarch64 Linux) for the ARMPL variant
## of the portable bundle. ARMPL is the BLAS/LAPACK the upstream Arm builds
## use, and the sources have dedicated code paths for it (INLA_WITH_ARMPL).
##
## Knobs:
##   ARMPL_VERSION   default below
##   ARMPL_PKG       rpm (RHEL/Alma/manylinux) or deb (Debian/Ubuntu)
##   ARMPL_URL       full tarball URL, when the permalink pattern changes
set -e

ARMPL_VERSION=${ARMPL_VERSION:-26.07}
ARMPL_PKG=${ARMPL_PKG:-rpm}
ARMPL_URL=${ARMPL_URL:-https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_${ARMPL_VERSION}/arm-performance-libraries_${ARMPL_VERSION}_${ARMPL_PKG}_gcc.tar}

if ls -d /opt/arm/armpl_* >/dev/null 2>&1; then
    echo "ARMPL already installed: $(ls -d /opt/arm/armpl_* | tail -1)"
    exit 0
fi

echo "fetching $ARMPL_URL"
curl -fL -o /tmp/armpl.tar "$ARMPL_URL"
mkdir -p /tmp/armpl && tar xf /tmp/armpl.tar -C /tmp/armpl --strip-components=1

## Two shapes have shipped over the years: a self-extracting installer
## script (-a accepts the licence non-interactively, -f overwrites, -i
## chooses the prefix), or plain packages to install directly.
INSTALLER=$(ls /tmp/armpl/*.sh 2>/dev/null | head -1 || true)
if [ -n "$INSTALLER" ]; then
    chmod +x "$INSTALLER"
    "$INSTALLER" -a -f -i /opt/arm
elif ls /tmp/armpl/*.rpm >/dev/null 2>&1; then
    rpm -Uvh --nodeps /tmp/armpl/*.rpm
elif ls /tmp/armpl/*.deb >/dev/null 2>&1; then
    dpkg -i /tmp/armpl/*.deb
else
    echo "ERROR: unrecognised ARMPL tarball layout:"; ls -R /tmp/armpl | head -20
    exit 1
fi

ARMPL_DIR=$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)
[ -n "$ARMPL_DIR" ] || { echo "ERROR: ARMPL installer produced no /opt/arm/armpl_* tree"; exit 1; }
ls "$ARMPL_DIR"/lib/libarmpl* >/dev/null || { echo "ERROR: no libarmpl in $ARMPL_DIR"; exit 1; }
echo "OK: ARMPL installed at $ARMPL_DIR"
