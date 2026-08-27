# quad-junctions (vendored subset)

This directory is a **trimmed, vendored copy** of the quad-junctions project
(https://github.com/TianyueLi75/quad-junctions), reduced to only the code that
GlymphBIE's capped-annulus geometry actually needs. It is no longer a git
submodule — the files are checked into this repo directly.

## What is here

- `include/sctl/experimental/` — the developed SCTL quadrilateral-element code
  (`quad_element.{hpp,cpp}`, `alpert_quadr.cpp`, `bench_quad.hpp`). This
  intentionally **shadows** the older ancestor shipped inside upstream SCTL, so
  the include order in the top-level `CMakeLists.txt` (`QJ_INCLUDE_DIR` before
  `SCTL_INCLUDE_DIR`) is load-bearing.
- `include/quad_junctions/`
  - `collar_mount_geom.hpp` — the butterfly / O-grid tip-cap mesher
    (`add_tip_cap_butterfly`, `flip_group`, `report_area`) used to close each
    slender tube into a watertight capsule.
  - `quad_scheme.hpp` — `QJDefaultScheme` / near-quadrature scheme selection for
    `QuadElemList::SetQuadScheme`.
  - `csbq_sctl_compat.hpp` — force-included (`-include`) CSBQ/SCTL shim
    (`SCTL_QUOTEME`).
  - `mpi_utils.hpp` — thin MPI helpers pulled in by `collar_mount_geom.hpp`.

## What was removed

The full quad-junctions repo (its `src/`, `python/`, `scripts/`, `bench/`, its
other `quad_junctions/*.hpp` geometry generators — cilia carpets, Y-bifurcations,
vessel networks, etc., `include/stokes_bio.*`, and its own pvfmm/CSBQ
submodules) is **not** vendored. None of it is referenced by GlymphBIE. See the
capped-annulus consumers: `include/CappedAnnulusGeom.hpp`,
`include/CappedAnnulusFlow.hpp`, and `examples/CappedAnnulus_*.cpp`.

To update this subset, pull the corresponding files from the upstream repo and
re-run the capped-annulus build + tests (`CappedAnnulus_watertight`,
`CappedAnnulus_SpatialSolver_tests`).
