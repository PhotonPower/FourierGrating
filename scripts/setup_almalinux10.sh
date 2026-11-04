#!/usr/bin/env bash
# Installs build dependencies of FourierGrating on AlmaLinux 10 (run as root).
#   ./scripts/setup_almalinux10.sh              CPU toolchain, libraries, Python, LaTeX
#   ./scripts/setup_almalinux10.sh --with-cuda  additionally the CUDA toolkit (NVIDIA repo rhel10)
# Missing optional packages are skipped; CMake falls back to FetchContent for Eigen/JSON/GTest.
set -euo pipefail

WITH_CUDA=0
for arg in "$@"; do
  case "$arg" in
    --with-cuda) WITH_CUDA=1 ;;
    -h|--help) sed -n '2,6p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 1 ;;
  esac
done

if [[ $EUID -ne 0 ]]; then echo "please run as root (sudo)" >&2; exit 1; fi

dnf -y install dnf-plugins-core epel-release
dnf config-manager --set-enabled crb || true
dnf -y makecache

dnf -y install --skip-unavailable \
  gcc-c++ make cmake ninja-build git pkgconf-pkg-config \
  eigen3-devel json-devel gtest-devel gmock-devel \
  flexiblas-devel lapack-devel fftw-devel hdf5-devel \
  python3 python3-numpy python3-matplotlib \
  texlive-scheme-basic texlive-collection-latexrecommended texlive-collection-langgerman \
  texlive-tcolorbox texlive-enumitem texlive-mathtools texlive-lm

if [[ $WITH_CUDA -eq 1 ]]; then
  dnf config-manager --add-repo \
    https://developer.download.nvidia.com/compute/cuda/repos/rhel10/x86_64/cuda-rhel10.repo
  dnf -y makecache
  dnf -y install cuda-toolkit
  cat << 'MSG'

CUDA toolkit installed (usually /usr/local/cuda). Add to your shell profile:
  export PATH=/usr/local/cuda/bin:$PATH
  export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
The NVIDIA driver is NOT installed by this script; follow the NVIDIA driver installation
guide for RHEL/AlmaLinux 10 and verify with `nvidia-smi`. Required: CUDA >= 12.6 Update 2.
MSG
fi

echo "done. Next: cmake --preset cpu-release && cmake --build --preset cpu-release && ctest --preset cpu-release"
