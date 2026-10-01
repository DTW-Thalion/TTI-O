#!/usr/bin/env bash
# setup-objc-wsl.sh — install the ObjC build environment on Ubuntu / WSL.
#
# Mirrors the CI composite actions (.github/actions/setup-gnustep,
# setup-hdf5, setup-libarrow) so a local machine builds and tests the
# ObjC SDK the way the objc-build-test job does: clang + libobjc2 +
# gnustep-make + gnustep-base from source, HDF5 from source, and
# libarrow-dev from Apache's apt source. Installs under /usr/local and
# needs sudo.
#
# Usage:
#   scripts/setup-objc-wsl.sh
#
# Afterwards:
#   . /usr/local/share/GNUstep/Makefiles/GNUstep.sh
#   cd objc && ./build.sh check

set -euo pipefail

LIBOBJC2_REF=${LIBOBJC2_REF:-v2.3}
TOOLS_MAKE_REF=${TOOLS_MAKE_REF:-make-2_9_3}
LIBS_BASE_REF=${LIBS_BASE_REF:-base-1_31_1}
HDF5_VERSION=${HDF5_VERSION:-1.14.6}
PREFIX=/usr/local
REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  clang cmake ninja-build make git ca-certificates zlib1g-dev libssl-dev \
  libxml2-dev libgnutls28-dev libffi-dev libicu-dev libblocksruntime-dev \
  libcurl4-openssl-dev libwebsockets-dev tzdata libzstd-dev samtools \
  lsb-release wget

if [ ! -e "$PREFIX/lib/libobjc.so" ]; then
  git clone --depth 1 --branch "$LIBOBJC2_REF" https://github.com/gnustep/libobjc2.git "$WORK/libobjc2"
  (cd "$WORK/libobjc2" && git submodule update --init --recursive \
    && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
    && ninja -C build && sudo ninja -C build install)
fi

if [ ! -e "$PREFIX/share/GNUstep/Makefiles/GNUstep.sh" ]; then
  git clone --depth 1 --branch "$TOOLS_MAKE_REF" https://github.com/gnustep/tools-make.git "$WORK/tools-make"
  (cd "$WORK/tools-make" && ./configure --prefix="$PREFIX" --with-library-combo=ng-gnu-gnu CC=clang \
    && make && sudo make install)
fi

if ! ls "$PREFIX"/lib/libgnustep-base.so* >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  . "$PREFIX/share/GNUstep/Makefiles/GNUstep.sh"
  git clone --depth 1 --branch "$LIBS_BASE_REF" https://github.com/gnustep/libs-base.git "$WORK/libs-base"
  (cd "$WORK/libs-base" && ./configure CC=clang && make && sudo -E make install)
fi
sudo ldconfig

if ! ls "$PREFIX"/lib/libhdf5.so* >/dev/null 2>&1; then
  bash "$REPO/scripts/install-hdf5.sh" "$HDF5_VERSION"
  sudo ldconfig "$PREFIX/lib"
fi

if ! pkg-config --exists arrow 2>/dev/null; then
  codename=$(lsb_release --codename --short)
  if wget -qO "$WORK/arrow-apt-source.deb" \
       "https://apache.jfrog.io/artifactory/arrow/ubuntu/apache-arrow-apt-source-latest-${codename}.deb"; then
    sudo apt-get install -y --no-install-recommends "$WORK/arrow-apt-source.deb"
    sudo apt-get update -qq
    sudo apt-get install -y --no-install-recommends libarrow-dev
  else
    echo "warning: no Apache Arrow apt source for '${codename}'; the transport" >&2
    echo "         tabular packets (TTIOArrowIpcBridge) will not build." >&2
  fi
fi

echo
echo "Done. In each new shell:"
echo "  . $PREFIX/share/GNUstep/Makefiles/GNUstep.sh"
echo "  cd objc && ./build.sh check"
