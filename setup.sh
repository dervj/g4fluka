#!/usr/bin/env bash
# setup.sh — configure the build environment and (re)compile g4fluka.
#
# Source this file (don't execute it) so that the env vars persist:
#   source setup.sh
# Or just run it to build:
#   bash setup.sh

set -euo pipefail

G4INSTALL=/cvmfs/muoncollider.cern.ch/release/2.9/linux-almalinux9-x86_64/gcc-11.4.1/geant4-11.2.0-gb5ve6rbe4ru5cglgtbnzh3xhjeejczd
G4CMAKE=${G4INSTALL}/lib64/cmake/Geant4

if [[ ! -d "$G4INSTALL" ]]; then
    echo "ERROR: Geant4 not found at $G4INSTALL"
    echo "       Is the muoncollider CVMFS repo mounted?"
    exit 1
fi

# Set all dataset environment variables (G4NEUTRONHPDATA, G4LEDATA, …)
source "${G4INSTALL}/bin/geant4.sh"

# Add Geant4 libraries to the linker path
export LD_LIBRARY_PATH="${G4INSTALL}/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "Geant4 $(${G4INSTALL}/bin/geant4-config --version) environment loaded."

# ---- Build --------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake .. -DGeant4_DIR="$G4CMAKE"
make -j"$(nproc)"

echo ""
echo "Build complete.  Executable: ${BUILD_DIR}/g4fluka"
echo ""
echo "Energy scans:"
echo "  cd ${BUILD_DIR} && ./g4fluka scan_tungsten.mac"
echo "  cd ${BUILD_DIR} && ./g4fluka scan_BPE.mac"
