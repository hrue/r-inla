#!/usr/bin/env bash
## Install Arm Performance Libraries (aarch64 Linux) for the ARMPL variant
## of the portable bundle. ARMPL is the BLAS/LAPACK the upstream Arm builds
## use, and the sources have dedicated code paths for it (INLA_WITH_ARMPL).
##
## Knobs:
##   ARMPL_VERSION   default below
##   ARMPL_URL       full tarball URL, when the permalink pattern changes
##   ARMPL_DISTRO    installer flavour matching the container (RHEL-8 for
##                   manylinux_2_28, Ubuntu-22.04 on a plain runner, ...)
set -e

ARMPL_VERSION=${ARMPL_VERSION:-24.10}
ARMPL_DISTRO=${ARMPL_DISTRO:-RHEL-8}
ARMPL_URL=${ARMPL_URL:-https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_${ARMPL_VERSION}/arm-performance-libraries_${ARMPL_VERSION}_${ARMPL_DISTRO}_gcc.tar}

if ls -d /opt/arm/armpl_* >/dev/null 2>&1; then
    echo "ARMPL already installed: $(ls -d /opt/arm/armpl_* | tail -1)"
    exit 0
fi

echo "fetching $ARMPL_URL"
curl -fL -o /tmp/armpl.tar "$ARMPL_URL"
mkdir -p /tmp/armpl && tar xf /tmp/armpl.tar -C /tmp/armpl --strip-components=1

## The tarball ships a self-extracting installer; -a accepts the licence
## non-interactively, -f overwrites a previous install.
INSTALLER=$(ls /tmp/armpl/*.sh | head -1)
chmod +x "$INSTALLER"
"$INSTALLER" -a -f -i /opt/arm

ARMPL_DIR=$(ls -d /opt/arm/armpl_* 2>/dev/null | sort -V | tail -1 || true)
[ -n "$ARMPL_DIR" ] || { echo "ERROR: ARMPL installer produced no /opt/arm/armpl_* tree"; exit 1; }
ls "$ARMPL_DIR"/lib/libarmpl* >/dev/null || { echo "ERROR: no libarmpl in $ARMPL_DIR"; exit 1; }
echo "OK: ARMPL installed at $ARMPL_DIR"
