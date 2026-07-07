#!/usr/bin/env bash
# run_g4fluka.sh — HTCondor job wrapper. Runs on a worker node that has
# the muoncollider.cern.ch CVMFS repo mounted (see requirements in
# g4fluka.sub). Sources the same Geant4 install used to build the binary
# on ap23 — the absolute CVMFS path is identical on every site that mounts
# the repo, so the prebuilt executable runs without rebuilding.
set -euo pipefail

MACRO="$1"

G4INSTALL=/cvmfs/muoncollider.cern.ch/release/2.9/linux-almalinux9-x86_64/gcc-11.4.1/geant4-11.2.0-gb5ve6rbe4ru5cglgtbnzh3xhjeejczd
source "${G4INSTALL}/bin/geant4.sh"
export LD_LIBRARY_PATH="${G4INSTALL}/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

chmod +x g4fluka
./g4fluka "${MACRO}"
