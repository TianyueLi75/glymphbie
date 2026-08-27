#!/usr/bin/env bash
#
# run_tracer_particles.sh
#
# Serial driver for the wall-function-generated capped-annulus flow with tracer
# particles. Runs at a small, verification scale and saves the whole run as just
# TWO files in vis/: a single VTKHDF (.vtkhdf) time series for the moving wall
# geometry and one for the tracer cloud -- open either directly in ParaView.
#
# The capped path is fully serial (no PVFMM / no MPI). The VTKHDF writer needs
# the HDF5 C library, so load a serial hdf5 module before configuring/running:
#     module load hdf5
# If you have not configured the serial build yet:
#     cmake -B build-serial -DENABLE_PVFMM=OFF -DENABLE_MPI=OFF
#     cmake --build build-serial -j
#
# Args passed through to the binary:
#     [Nelem fourier capOrder capNaz Nsteps A Q tol]
# Defaults below are small (verification only). A is the peristaltic breathing
# amplitude; Q is the superimposed inflow/outflow (pin-pout) throughput flux.
set -euo pipefail

cd "$(dirname "$0")"

BIN=build-serial/bin/CappedAnnulus_TracerParticles_ex
if [[ ! -x "$BIN" ]]; then
  echo "error: $BIN not found. Configure and build the serial target first:" >&2
  echo "  cmake -B build-serial -DENABLE_PVFMM=OFF -DENABLE_MPI=OFF" >&2
  echo "  cmake --build build-serial -j" >&2
  exit 1
fi

mkdir -p vis

# Nelem fourier capOrder capNaz Nsteps A     Q     tol
NELEM=${1:-6}
FOURIER=${2:-16}
CAPORDER=${3:-8}
CAPNAZ=${4:-16}
NSTEPS=${5:-4}
A=${6:-0.02}
Q=${7:-0.001}
TOL=${8:-1e-8}

echo "Running $BIN  Nelem=$NELEM fourier=$FOURIER capOrder=$CAPORDER capNaz=$CAPNAZ Nsteps=$NSTEPS A=$A Q=$Q tol=$TOL"
"$BIN" "$NELEM" "$FOURIER" "$CAPORDER" "$CAPNAZ" "$NSTEPS" "$A" "$Q" "$TOL"

echo
echo "Output written to vis/ (open these .vtkhdf files directly in ParaView):"
ls -1 vis/walls.vtkhdf vis/tracers.vtkhdf 2>/dev/null || true
