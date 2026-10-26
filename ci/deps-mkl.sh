#!/usr/bin/env bash
## Install Intel oneAPI MKL, the BLAS/LAPACK the upstream Linux x86_64
## builds link statically. Works on both apt (Ubuntu runner) and dnf/yum
## (the manylinux container) hosts.
set -e

SUDO=""
[ "$(id -u)" != 0 ] && SUDO=sudo

if [ -d /opt/intel/oneapi/mkl/latest ]; then
    echo "MKL already installed"
    exit 0
fi

if command -v apt-get >/dev/null 2>&1; then
    curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
        | gpg --dearmor | $SUDO tee /usr/share/keyrings/oneapi.gpg > /dev/null
    echo "deb [signed-by=/usr/share/keyrings/oneapi.gpg] https://apt.repos.intel.com/oneapi all main" \
        | $SUDO tee /etc/apt/sources.list.d/oneAPI.list > /dev/null
    $SUDO apt-get update
    $SUDO apt-get install -y intel-oneapi-mkl-devel
else
    $SUDO tee /etc/yum.repos.d/oneAPI.repo > /dev/null <<'EOF'
[oneAPI]
name=Intel oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
    $SUDO dnf -y install intel-oneapi-mkl-devel
fi

ls -d /opt/intel/oneapi/mkl/latest >/dev/null || { echo "ERROR: MKL not installed"; exit 1; }
echo "OK: MKL at /opt/intel/oneapi/mkl/latest"
