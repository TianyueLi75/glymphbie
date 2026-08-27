#ifndef _CAPPEDANNULUSGEOM_HPP_
#define _CAPPEDANNULUSGEOM_HPP_

/**
 * CappedAnnulusGeom.hpp
 *
 * Shared geometry builders for the non-periodic, *capped* framing of the
 * annular perivascular channel. A single straight wall (inner arteriole or
 * outer glial wall) is a CSBQ SlenderElemList tube whose two axial ends are
 * closed by hemispherical butterfly / O-grid dome caps from quad-junctions
 * (add_tip_cap_butterfly), producing a watertight capsule. Two nested capsules
 * bound the perivascular fluid.
 *
 * These builders were originally written inline in
 * examples/CappedAnnulus_watertight.cpp; they are factored here so the
 * watertightness test and the flow drivers (accuracy / self-convergence /
 * wall-coupled time stepping) share one definition.
 *
 * Header-only; template on Real. Assumes a SCTL that provides
 * sctl::SlenderElemList and sctl::QuadElemList, plus the quad-junctions cap
 * mesher <quad_junctions/collar_mount_geom.hpp>. Include those before this file
 * (every capped driver does `#include <sctl.hpp>` and the quad_junctions
 * headers first).
 */

#include <sctl.hpp>
#include <quad_junctions/collar_mount_geom.hpp>   // add_tip_cap_butterfly, flip_group
#include <functional>
#include <utility>

namespace glymphbie {

// Degree-5 "smootherstep" partition-of-unity weight on [0,1]:
//   S(t) = t^3 (6 t^2 - 15 t + 10),   S(0)=0, S(1)=1, S'=S''=0 at both ends.
// This is the same quintic used by the quad-junctions transition blends
// (quad_junctions::pou_weight, HybridAssembly::pou_ramp, collar_s5); it is
// C^2 at the seams and, being degree 5, is represented EXACTLY by an order>=6
// Chebyshev panel -- so a CSBQ transition panel that samples it has no blend
// representation error. We keep a local copy rather than call
// quad_junctions::pou_weight to avoid depending on its mutable pou_kind() global.
template <class Real>
Real pou_s5(Real t) {
  if (t <= (Real)0) return (Real)0;
  if (t >= (Real)1) return (Real)1;
  return t * t * t * ((Real)6 * t * t - (Real)15 * t + (Real)10);
}

// Straight-tube (cylinder) body of revolution about the x-axis, centerline at
// (x, y0, z0) for x in [0, L], constant radius r. Ordered [x1,y1,z1, ...].
template <class Real>
sctl::SlenderElemList<Real> BuildTube(sctl::Long Nelem, sctl::Integer cheb, sctl::Integer fourier,
                                      Real L, Real r, Real y0, Real z0) {
  sctl::Vector<sctl::Long> elem_order, forder;
  sctl::Vector<Real> coord, radius;
  for (sctl::Long k = 0; k < Nelem; k++) {
    elem_order.PushBack(cheb);
    forder.PushBack(fourier);
    const sctl::Vector<Real>& cn = sctl::SlenderElemList<Real>::CenterlineNodes(cheb);
    for (sctl::Integer j = 0; j < cheb; j++) {
      const Real x = L * (k + cn[j]) / Nelem;
      coord.PushBack(x); coord.PushBack(y0); coord.PushBack(z0);
      radius.PushBack(r);
    }
  }
  return sctl::SlenderElemList<Real>(elem_order, forder, coord, radius);
}

// Straight-tube body of revolution about the x-axis (centerline at (x,y0,z0),
// x in [0,L]) whose radius follows an arbitrary axial profile r_of_x(x). Same
// node layout as BuildTube -- only the constant radius is replaced by the
// per-node profile value. Used to build the transition + moving-wall tube whose
// radius blends a fixed cap radius into a wall-function-defined radius.
template <class Real>
sctl::SlenderElemList<Real> BuildProfiledTube(sctl::Long Nelem, sctl::Integer cheb, sctl::Integer fourier,
                                              Real L, Real y0, Real z0,
                                              const std::function<Real(Real)>& r_of_x) {
  sctl::Vector<sctl::Long> elem_order, forder;
  sctl::Vector<Real> coord, radius;
  for (sctl::Long k = 0; k < Nelem; k++) {
    elem_order.PushBack(cheb);
    forder.PushBack(fourier);
    const sctl::Vector<Real>& cn = sctl::SlenderElemList<Real>::CenterlineNodes(cheb);
    for (sctl::Integer j = 0; j < cheb; j++) {
      const Real x = L * (k + cn[j]) / Nelem;
      coord.PushBack(x); coord.PushBack(y0); coord.PushBack(z0);
      radius.PushBack(r_of_x(x));
    }
  }
  return sctl::SlenderElemList<Real>(elem_order, forder, coord, radius);
}

// Blended axial radius profile r(x) and its time-derivative rdot(x) for a
// capped, moving-wall tube. Over the two transition bands [0, xt_lo] and
// [L - xt_hi, L] the radius blends the FIXED cap radius R_cap (constant in
// time) into the wall-function radius r_wall(x) via the degree-5 POU weight w;
// in the interior [xt_lo, L - xt_hi] it is exactly r_wall(x):
//   r(x)    = (1 - w) * R_cap + w * r_wall(x)
//   rdot(x) = w * rdot_wall(x)          (R_cap is time-constant)
// with w = pou_s5(eta), eta the in-band fraction (0 at the cap seam, 1 at the
// wall seam). At the cap seam w = w' = 0 so r = R_cap with zero slope -> the
// cap ring is fixed and watertight; at the wall seam w = 1, w' = 0 so r and its
// slope equal the wall's -> tangent-matched to the moving wall (this requires
// r_wall(x) to be evaluable across the transition band, i.e. as a smooth
// function of x, not just at the seam node).
template <class Real>
struct CappedRadiusProfile {
  std::function<Real(Real)> r;      // radius r(x)
  std::function<Real(Real)> rdot;   // radial velocity rdot(x)
};

template <class Real>
CappedRadiusProfile<Real> MakeCappedRadiusProfile(Real R_cap, Real xt_lo, Real xt_hi, Real L,
                                                  std::function<Real(Real)> r_wall,
                                                  std::function<Real(Real)> rdot_wall) {
  CappedRadiusProfile<Real> prof;
  // w(x): POU weight, 0 on the fixed cap, 1 on the moving wall.
  auto weight = [xt_lo, xt_hi, L](Real x) -> Real {
    if (x < xt_lo)      return pou_s5<Real>(x / xt_lo);          // low-x band: 0 at x=0
    if (x > L - xt_hi)  return pou_s5<Real>((L - x) / xt_hi);    // high-x band: 0 at x=L
    return (Real)1;                                             // interior wall region
  };
  prof.r = [R_cap, r_wall, weight](Real x) -> Real {
    const Real w = weight(x);
    return ((Real)1 - w) * R_cap + w * r_wall(x);
  };
  prof.rdot = [rdot_wall, weight](Real x) -> Real {
    return weight(x) * rdot_wall(x);
  };
  return prof;
}

// Orient a just-built cap group so its normals point radially OUTWARD from the
// tube (away from the tip point Ctip, i.e. into the fluid / out of the solid,
// matching the CSBQ shaft's un-flippable radially-outward normal). Accumulates
// n.(x - Ctip) over the cap's nodes and transposes each element (flip_group) if
// that is negative. Returns true if flipped.
template <class Real>
bool OrientCapOutward(sctl::Vector<Real>& Xcap, sctl::Integer order, const Real Ctip[3]) {
  sctl::QuadElemList<Real> tmp(order, Xcap);
  sctl::Vector<Real> Xc, Xnc; tmp.GetNodeCoord(&Xc, &Xnc, nullptr);
  Real acc = 0;
  for (sctl::Long i = 0; i < Xc.Dim() / 3; i++)
    for (int k = 0; k < 3; k++) acc += Xnc[i*3+k] * (Xc[i*3+k] - Ctip[k]);
  if (acc >= 0) return false;
  // Transpose each order x order element (swap u<->v) to negate its normals (cf. flip_group).
  const sctl::Long nn = (sctl::Long)order * order, ne = Xcap.Dim() / (nn * 3);
  for (sctl::Long e = 0; e < ne; e++)
    for (sctl::Integer i = 0; i < order; i++)
      for (sctl::Integer j = i + 1; j < order; j++)
        for (int c = 0; c < 3; c++)
          std::swap(Xcap[(e*nn + i*order + j)*3 + c], Xcap[(e*nn + j*order + i)*3 + c]);
  return true;
}

// Build both dome caps (at x=0 and x=L) for a tube of radius r about the x-axis
// through (*, y0, z0), oriented outward, as a single QuadElemList.
template <class Real>
sctl::QuadElemList<Real> BuildCaps(sctl::Integer order, sctl::Integer Naz, Real L, Real r,
                                   Real y0, Real z0, const sctl::Comm& comm) {
  const Real w1[3] = {0, 1, 0}, w2[3] = {0, 0, 1};  // frame perpendicular to the x-axis travel tangent
  sctl::Vector<Real> Xcaps;

  // Low-x cap: tip at (0, y0, z0), dome bulges toward -x.
  {
    const Real Ctip[3] = {(Real)0, y0, z0}, Ttip[3] = {(Real)-1, (Real)0, (Real)0};
    sctl::Vector<Real> Xc;
    quad_junctions::add_tip_cap_butterfly<Real>(Xc, order, Ctip, Ttip, w1, w2, r, Naz);
    OrientCapOutward<Real>(Xc, order, Ctip);
    for (sctl::Long i = 0; i < Xc.Dim(); i++) Xcaps.PushBack(Xc[i]);
  }
  // High-x cap: tip at (L, y0, z0), dome bulges toward +x.
  {
    const Real Ctip[3] = {L, y0, z0}, Ttip[3] = {(Real)1, (Real)0, (Real)0};
    sctl::Vector<Real> Xc;
    quad_junctions::add_tip_cap_butterfly<Real>(Xc, order, Ctip, Ttip, w1, w2, r, Naz);
    OrientCapOutward<Real>(Xc, order, Ctip);
    for (sctl::Long i = 0; i < Xc.Dim(); i++) Xcaps.PushBack(Xc[i]);
  }
  return sctl::QuadElemList<Real>(order, Xcaps, comm);
}

}  // namespace glymphbie

#endif  // _CAPPEDANNULUSGEOM_HPP_
