#ifndef _FLOWVIS_HPP_
#define _FLOWVIS_HPP_

/**
 * FlowVis.hpp
 *
 * Small helpers for visualizing the annular perivascular flow field as a VTK
 * time series, alongside the Lagrangian tracers of Particle.hpp.
 *
 * A flow-field snapshot is just a set of interior sample points colored by the
 * evaluated velocity -- structurally identical to a tracer point cloud -- so we
 * reuse TracerParticles<Real>::WriteVTK to render it (no new VTU code):
 *   - MakeAnnularGrid builds a small cylindrical grid of points strictly inside
 *     the annular fluid gap (filtered by TracerParticles::InDomain, so no target
 *     lands on/near a surface where the BIE evaluation is near-singular),
 *   - WriteFieldVTK wraps those points + a velocity field in a throwaway
 *     TracerParticles and writes them as one VTU point cloud.
 *
 * Header-only; template on Real. The velocity itself is produced by the caller
 * (e.g. CappedAnnulusFlow::EvaluateVelocity) -- this header is solver-agnostic.
 */

#include <sctl.hpp>
#include "Particle.hpp"
#include <string>
#include <cmath>

namespace glymphbie {

// Build a small cylindrical grid of sample points strictly inside the annular
// fluid gap of geometry g: Nx axial stations in (x0,x1), Nr radial shells in
// (R_in,R_out), Nth azimuthal angles. Points that fail the in-domain test (the
// same radial/axial band test the tracers use, with g.in_tol margin) are
// dropped, so the returned AoS coords {x1,y1,z1,...} are safe off-surface
// targets for a BIE velocity evaluation.
template <class Real>
sctl::Vector<Real> MakeAnnularGrid(const typename TracerParticles<Real>::Geometry& g,
                                   sctl::Long Nx, sctl::Long Nr, sctl::Long Nth) {
  const Real L = g.x1 - g.x0;
  const Real pi = sctl::const_pi<Real>();
  sctl::Vector<Real> X;
  for (sctl::Long ix = 0; ix < Nx; ix++) {
    const Real x = g.x0 + L * (Real)(ix + 1) / (Real)(Nx + 1);          // interior axial stations
    for (sctl::Long ir = 0; ir < Nr; ir++) {
      const Real r = g.R_in + (g.R_out - g.R_in) * (Real)(ir + 1) / (Real)(Nr + 1);  // interior shells
      for (sctl::Long it = 0; it < Nth; it++) {
        const Real th = 2 * pi * (Real)it / (Real)Nth;
        X.PushBack(x);
        X.PushBack(g.y0 + r * std::cos(th));
        X.PushBack(g.z0 + r * std::sin(th));
      }
    }
  }
  // Keep only strictly-interior points (defensive: the shells above are already
  // interior, but this guards against margin/round-off and any future layout).
  TracerParticles<Real> probe(X, g);
  const sctl::Vector<sctl::Long> in = probe.InDomain();
  sctl::Vector<Real> Xin;
  for (sctl::Long i = 0; i < in.Dim(); i++)
    if (in[i]) for (sctl::Integer k = 0; k < 3; k++) Xin.PushBack(X[i*3+k]);
  return Xin;
}

// Write a flow-field snapshot: the points X (AoS) colored by the velocity U
// (AoS, same length), as one VTU point cloud. `fname` gets ".vtu". Reuses the
// tracer point-cloud writer.
template <class Real>
void WriteFieldVTK(const std::string& fname, const sctl::Vector<Real>& X,
                   const sctl::Vector<Real>& U,
                   const typename TracerParticles<Real>::Geometry& g = typename TracerParticles<Real>::Geometry()) {
  TracerParticles<Real> field(X, g);
  field.WriteVTK(fname, U);
}

}  // namespace glymphbie

#endif  // _FLOWVIS_HPP_
