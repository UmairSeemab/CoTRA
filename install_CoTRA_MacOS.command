#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

echo "============================================================"
echo " CoTRA native installer for macOS"
echo " Intel and Apple Silicon (M1/M2/M3/M4...)"
echo "============================================================"
echo

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript was not found."
  echo "Install R first, then run this installer again."
  echo "https://cran.r-project.org/bin/macosx/"
  echo
  read -r -p "Press Enter to close..." || true
  exit 1
fi

if ! xcode-select -p >/dev/null 2>&1; then
  echo "Apple Command Line Tools are required."
  echo "Starting Apple's installer..."
  xcode-select --install || true
  echo
  echo "After the Command Line Tools installation finishes, run this file again."
  read -r -p "Press Enter to close..." || true
  exit 1
fi

# GUI-launched applications do not always inherit Homebrew's PATH.
if [[ -x /opt/homebrew/bin/brew ]]; then
  BREW=/opt/homebrew/bin/brew
elif [[ -x /usr/local/bin/brew ]]; then
  BREW=/usr/local/bin/brew
elif command -v brew >/dev/null 2>&1; then
  BREW="$(command -v brew)"
else
  echo "ERROR: Homebrew was not found."
  echo
  echo "Install Homebrew from:"
  echo "https://brew.sh"
  echo
  echo "Then run this installer again."
  echo "CoTRA uses Homebrew HDF5, pkgconf and OpenSSL when R packages"
  echo "must be compiled from source on macOS."
  read -r -p "Press Enter to close..." || true
  exit 1
fi

echo "Using Homebrew: $BREW"
echo "Installing/checking required macOS libraries..."
"$BREW" install hdf5 pkgconf openssl@3

BREW_PREFIX="$("$BREW" --prefix)"
HDF5_PREFIX="$("$BREW" --prefix hdf5)"
OPENSSL_PREFIX="$("$BREW" --prefix openssl@3)"
PKGCONF_PREFIX="$("$BREW" --prefix pkgconf)"

export PATH="$HDF5_PREFIX/bin:$PKGCONF_PREFIX/bin:$BREW_PREFIX/bin:$PATH"
export PKG_CONFIG_PATH="$HDF5_PREFIX/lib/pkgconfig:$OPENSSL_PREFIX/lib/pkgconfig:$PKGCONF_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPATH="$OPENSSL_PREFIX/include:$HDF5_PREFIX/include:${CPATH:-}"
export LIBRARY_PATH="$OPENSSL_PREFIX/lib:$HDF5_PREFIX/lib:${LIBRARY_PATH:-}"
export CPPFLAGS="-I$OPENSSL_PREFIX/include -I$HDF5_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$OPENSSL_PREFIX/lib -L$HDF5_PREFIX/lib ${LDFLAGS:-}"

echo
echo "macOS prerequisites detected:"
echo "  Architecture: $(uname -m)"
echo "  HDF5:        $HDF5_PREFIX"
echo "  OpenSSL:     $OPENSSL_PREFIX"
echo "  pkg-config:  $(command -v pkg-config || true)"
echo "  h5cc:        $(command -v h5cc || true)"
echo

if ! command -v pkg-config >/dev/null 2>&1; then
  echo "ERROR: pkg-config is still unavailable after installing pkgconf."
  read -r -p "Press Enter to close..." || true
  exit 1
fi

if ! command -v h5cc >/dev/null 2>&1; then
  echo "ERROR: h5cc is still unavailable after installing HDF5."
  read -r -p "Press Enter to close..." || true
  exit 1
fi

Rscript install_cotra_packages.R

echo
echo "CoTRA installation finished."
echo "Open R/RStudio and run:"
echo "  library(CoTRA)"
echo "  CoTRA::runCoTRA()"
echo
read -r -p "Press Enter to close..." || true
