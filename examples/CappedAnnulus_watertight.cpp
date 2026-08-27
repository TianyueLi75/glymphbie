/**
 * CappedAnnulus_watertight.cpp
 *
 * Parallel-in-spirit companion to SpatialSolver_tests.cpp. Where that driver
 * closes the annular perivascular domain by *periodicity* at the two ends of
 * the channel, this one closes it *geometrically*: each of the two slender
 * walls (the inner arteriole wall and the outer glial wall -- each a CSBQ
 * SlenderElemList "arm") is capped at both ends by a butterfly / O-grid
 * hemispherical dome from quad-junctions (add_tip_cap_butterfly). The result is
 * two nested, watertight capsules; the perivascular fluid is the closed cavity
 * between them.
 *
 * There is no combined annulus mesh: the caps are mounted on the inner and
 * outer walls *separately* (they are ordinary bodies of revolution, so the
 * quad-junctions tip-cap mount applies directly to each).
 *
 * The program is a *geometry / watertightness* test -- no BIE solve, no PVFMM,
 * no MPI required. For each capsule it accumulates, over the far-field
 * quadrature of the tube (SlenderElemList) + its two caps (QuadElemList):
 *     area   = sum(wts)                             (compare 2*pi*r*L + 4*pi*r^2)
 *     flux   = sum(wts * n)      -> ~0 for closed   (divergence of a constant)
 *     volume = sum(wts * x.n)/3  = enclosed volume  (divergence theorem)
 * matching quad-junctions' report_area / areatest.cpp closure checks. A
 * watertight capsule has |flux|/area ~ 0 and volume ~ pi*r^2*L + (4/3)*pi*r^3.
 *
 * Build (serial, no PVFMM/MPI):
 *   see CMakeLists.txt target `CappedAnnulus_watertight` (ENABLE_PVFMM=OFF,
 *   ENABLE_MPI=OFF is enough).
 *
 * Usage: ./CappedAnnulus_watertight [Nelem ElemOrder FourierOrder CapNaz]
 */
#include <sctl.hpp>
#include <quad_junctions/collar_mount_geom.hpp>   // add_tip_cap_butterfly, flip_group, report_area
#include <quad_junctions/quad_scheme.hpp>         // QJDefaultScheme
#include <csbq/slender_element.hpp>
#include <csbq/slender_element.cpp>
#include <CappedAnnulusGeom.hpp>   // shared pou_s5 / BuildProfiledTube / MakeCappedRadiusProfile
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace sctl;
using namespace quad_junctions;

// Straight-tube (cylinder) body of revolution about the x-axis, centerline at
// (x, y0, z0) for x in [0, L], constant radius r. Ordered [x1,y1,z1, ...].
template <class Real>
static SlenderElemList<Real> BuildTube(Long Nelem, Integer cheb, Integer fourier, Real L, Real r, Real y0, Real z0) {
  Vector<Long> elem_order, forder;
  Vector<Real> coord, radius;
  for (Long k = 0; k < Nelem; k++) {
    elem_order.PushBack(cheb);
    forder.PushBack(fourier);
    const Vector<Real>& cn = SlenderElemList<Real>::CenterlineNodes(cheb);
    for (Integer j = 0; j < cheb; j++) {
      const Real x = L * (k + cn[j]) / Nelem;
      coord.PushBack(x); coord.PushBack(y0); coord.PushBack(z0);
      radius.PushBack(r);
    }
  }
  return SlenderElemList<Real>(elem_order, forder, coord, radius);
}

// Orient a just-built cap group so its normals point radially OUTWARD from the
// tube (away from the tip point Ctip, i.e. into the fluid / out of the solid,
// matching the CSBQ shaft's un-flippable radially-outward normal). Accumulates
// n.(x - Ctip) over the cap's nodes and transposes each element (flip_group) if
// that is negative. Returns true if flipped.
template <class Real>
static bool OrientCapOutward(Vector<Real>& Xcap, Integer order, const Real Ctip[3]) {
  QuadElemList<Real> tmp(order, Xcap);
  Vector<Real> Xc, Xnc; tmp.GetNodeCoord(&Xc, &Xnc, nullptr);
  Real acc = 0;
  for (Long i = 0; i < Xc.Dim() / 3; i++)
    for (int k = 0; k < 3; k++) acc += Xnc[i*3+k] * (Xc[i*3+k] - Ctip[k]);
  if (acc >= 0) return false;
  // Transpose each order x order element (swap u<->v) to negate its normals (cf. flip_group).
  const Long nn = (Long)order * order, ne = Xcap.Dim() / (nn * 3);
  for (Long e = 0; e < ne; e++)
    for (Integer i = 0; i < order; i++)
      for (Integer j = i + 1; j < order; j++)
        for (int c = 0; c < 3; c++)
          std::swap(Xcap[(e*nn + i*order + j)*3 + c], Xcap[(e*nn + j*order + i)*3 + c]);
  return true;
}

// Build both dome caps (at x=0 and x=L) for a tube of radius r about the x-axis
// through (*, y0, z0), oriented outward, as a single QuadElemList.
template <class Real>
static QuadElemList<Real> BuildCaps(Integer order, Integer Naz, Real L, Real r, Real y0, Real z0, const Comm& comm) {
  const Real w1[3] = {0, 1, 0}, w2[3] = {0, 0, 1};  // frame perpendicular to the x-axis travel tangent
  Vector<Real> Xcaps;

  // Low-x cap: tip at (0, y0, z0), dome bulges toward -x.
  {
    const Real Ctip[3] = {(Real)0, y0, z0}, Ttip[3] = {(Real)-1, (Real)0, (Real)0};
    Vector<Real> Xc;
    add_tip_cap_butterfly<Real>(Xc, order, Ctip, Ttip, w1, w2, r, Naz);
    OrientCapOutward<Real>(Xc, order, Ctip);
    for (Long i = 0; i < Xc.Dim(); i++) Xcaps.PushBack(Xc[i]);
  }
  // High-x cap: tip at (L, y0, z0), dome bulges toward +x.
  {
    const Real Ctip[3] = {L, y0, z0}, Ttip[3] = {(Real)1, (Real)0, (Real)0};
    Vector<Real> Xc;
    add_tip_cap_butterfly<Real>(Xc, order, Ctip, Ttip, w1, w2, r, Naz);
    OrientCapOutward<Real>(Xc, order, Ctip);
    for (Long i = 0; i < Xc.Dim(); i++) Xcaps.PushBack(Xc[i]);
  }
  return QuadElemList<Real>(order, Xcaps, comm);
}

// Accumulate area / flux / (x.n) closure sums over one element list (works for
// both SlenderElemList and QuadElemList via GetFarFieldNodes -- pure quadrature,
// no FMM). Mirrors report_area / areatest.cpp accum.
template <class Real, class LST>
static void AccumClosure(const LST& lst, Real tol, Real& area, Real flux[3], Real& xdotn, Long& nnode) {
  Vector<Real> X, Xn, wts, dist_far; Vector<Long> cnt;
  lst.GetFarFieldNodes(X, Xn, wts, dist_far, cnt, tol);
  for (Long i = 0; i < wts.Dim(); i++) {
    area += wts[i];
    for (int k = 0; k < 3; k++) flux[k] += wts[i] * Xn[i*3+k];
    for (int k = 0; k < 3; k++) xdotn  += wts[i] * X[i*3+k] * Xn[i*3+k];
  }
  nnode += X.Dim() / 3;
}

// Report the watertightness of one capsule = tube + its two caps.
template <class Real>
static Real ReportCapsule(const char* tag, const SlenderElemList<Real>& tube, const QuadElemList<Real>& caps,
                          Real L, Real r, Real tol, const Comm& comm) {
  Real area = 0, flux[3] = {0, 0, 0}, xdotn = 0; Long nnode = 0;
  AccumClosure<Real>(tube, tol, area, flux, xdotn, nnode);
  AccumClosure<Real>(caps, tol, area, flux, xdotn, nnode);

  area  = GlobalReduce((double)area,  comm, CommOp::SUM);
  xdotn = GlobalReduce((double)xdotn, comm, CommOp::SUM);
  for (int k = 0; k < 3; k++) flux[k] = GlobalReduce((double)flux[k], comm, CommOp::SUM);
  nnode = GlobalReduce((Long)nnode, comm, CommOp::SUM);

  const Real pi = const_pi<Real>();
  const Real vol = xdotn / 3;
  const Real flux_mag = sqrt<Real>(flux[0]*flux[0] + flux[1]*flux[1] + flux[2]*flux[2]);
  const Real area_exact = 2*pi*r*L + 4*pi*r*r;              // cylinder wall + two hemispheres
  const Real vol_exact  = pi*r*r*L + (Real)(4.0/3.0)*pi*r*r*r;  // cylinder + sphere
  if (!comm.Rank()) {
    std::cout << std::setprecision(10);
    std::cout << "  [" << tag << "] nodes=" << nnode << "  radius=" << r << "  length=" << L << "\n"
              << "      surface area = " << area << "   (exact " << area_exact
              << ",  rel err " << std::setprecision(3) << std::fabs(area - area_exact)/area_exact << ")\n"
              << std::setprecision(10)
              << "      volume       = " << vol << "   (exact " << vol_exact
              << ",  rel err " << std::setprecision(3) << std::fabs(vol - vol_exact)/vol_exact << ")\n"
              << "      |int n dA|   = " << flux_mag << "   (rel " << flux_mag / area
              << ")   <- watertightness (want ~0)\n";
  }
  return vol;
}

// Watertightness of a MOVING-WALL (radius-profiled) capsule: a blended
// fixed-cap -> wall-function tube (glymphbie::BuildProfiledTube) closed by the
// FIXED dome caps at R_cap. Returns |int n dA| / area (want ~0). The exact
// area/volume are shape-dependent here, so we report only the shape-independent
// closure ratio -- a nonzero value would mean the transition->cap join opened.
template <class Real>
static Real ReportProfiledWatertight(const char* tag, Long Nelem, Integer cheb, Integer fourier,
                                     Integer capOrd, Integer capNaz, Real L, Real R_cap,
                                     Real y0, Real z0, const std::function<Real(Real)>& r_of_x,
                                     Real tol, const Comm& comm) {
  SlenderElemList<Real> tube = glymphbie::BuildProfiledTube<Real>(Nelem, cheb, fourier, L, y0, z0, r_of_x);
  QuadElemList<Real>    caps = glymphbie::BuildCaps<Real>(capOrd, capNaz, L, R_cap, y0, z0, comm);
  caps.SetQuadScheme(QJDefaultScheme<Real>(), 10, 200, 12);

  Real area = 0, flux[3] = {0, 0, 0}, xdotn = 0; Long nnode = 0;
  AccumClosure<Real>(tube, tol, area, flux, xdotn, nnode);
  AccumClosure<Real>(caps, tol, area, flux, xdotn, nnode);
  area = GlobalReduce((double)area, comm, CommOp::SUM);
  for (int k = 0; k < 3; k++) flux[k] = GlobalReduce((double)flux[k], comm, CommOp::SUM);
  const Real flux_mag = sqrt<Real>(flux[0]*flux[0] + flux[1]*flux[1] + flux[2]*flux[2]);
  const Real ratio = flux_mag / area;
  if (!comm.Rank())
    std::cout << "  [" << tag << "]  R_cap=" << R_cap << "  |int n dA|/area = "
              << std::scientific << std::setprecision(3) << ratio << std::defaultfloat
              << "   <- watertightness across wall motion (want ~0)\n";
  return ratio;
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  using Real = double;
  {
    const Comm comm = Comm::World();

    // Discretization (positional args mirror the other examples' Nelem ElemOrder FourierOrder).
    const Long    Nelem   = (argc > 1) ? std::atol(argv[1]) : 8;
    const Integer cheb    = (argc > 2) ? std::atoi(argv[2]) : 10;
    const Integer fourier = (argc > 3) ? std::atoi(argv[3]) : 16;
    const Integer capNaz  = (argc > 4) ? std::atoi(argv[4]) : 16;  // azimuthal panels around the cap equator
    const Integer capOrd  = cheb;                                  // reuse Chebyshev order for the cap patches
    const Real    tol     = 1e-10;

    // Annular geometry: inner + outer straight walls, centerline along x in [0,L] at (y,z)=(0.5,0.5).
    const Real L = 1.0, y0 = 0.5, z0 = 0.5;
    const Real r_inner = 0.10;   // arteriole wall
    const Real r_outer = 0.30;   // glial wall

    if (!comm.Rank())
      std::cout << "Capped-annulus watertightness test  (Nelem=" << Nelem << " cheb=" << cheb
                << " fourier=" << fourier << " capNaz=" << capNaz << ")\n"
                << "  inner + outer slender walls, each closed by two butterfly dome caps.\n";

    // Inner capsule.
    SlenderElemList<Real> inner_tube = BuildTube<Real>(Nelem, cheb, fourier, L, r_inner, y0, z0);
    QuadElemList<Real>    inner_caps = BuildCaps<Real>(capOrd, capNaz, L, r_inner, y0, z0, comm);
    inner_caps.SetQuadScheme(QJDefaultScheme<Real>(), 10, 200, 12);

    // Outer capsule.
    SlenderElemList<Real> outer_tube = BuildTube<Real>(Nelem, cheb, fourier, L, r_outer, y0, z0);
    QuadElemList<Real>    outer_caps = BuildCaps<Real>(capOrd, capNaz, L, r_outer, y0, z0, comm);
    outer_caps.SetQuadScheme(QJDefaultScheme<Real>(), 10, 200, 12);

    const Real vol_in  = ReportCapsule<Real>("inner capsule", inner_tube, inner_caps, L, r_inner, tol, comm);
    const Real vol_out = ReportCapsule<Real>("outer capsule", outer_tube, outer_caps, L, r_outer, tol, comm);

    if (!comm.Rank()) {
      std::cout << std::setprecision(10)
                << "  [perivascular cavity]  fluid volume = outer - inner = "
                << (vol_out - vol_in) << "\n";
    }

    // --- Moving-wall watertightness: fixed caps + blended transition tube ---
    // At a nonzero peristaltic wall state r_wall(x) = R_cap + A sin(2 pi x/L),
    // the transition blends R_cap (at the caps) into r_wall (interior). A
    // watertight |int n dA|/area ~ 0 confirms the transition->cap seam stays
    // sealed even though the interior wall has moved off the cap radius.
    if (!comm.Rank()) std::cout << "  --- moving-wall (profiled) watertightness ---\n";
    {
      // Regression guard: a CONSTANT profile r(x) = R must reproduce BuildTube(R)
      // node-for-node, so the moving-wall path degenerates to the static path.
      SlenderElemList<Real> tube_const = BuildTube<Real>(Nelem, cheb, fourier, L, r_inner, y0, z0);
      SlenderElemList<Real> tube_prof  = glymphbie::BuildProfiledTube<Real>(
          Nelem, cheb, fourier, L, y0, z0, [r_inner](Real){ return r_inner; });
      Vector<Real> Xa, Xb; tube_const.GetNodeCoord(&Xa, nullptr, nullptr);
      tube_prof.GetNodeCoord(&Xb, nullptr, nullptr);
      Real dmax = 0; for (Long i = 0; i < Xa.Dim(); i++) dmax = std::max(dmax, std::fabs(Xa[i]-Xb[i]));
      if (!comm.Rank())
        std::cout << "  [const-profile regression]  max |BuildProfiledTube - BuildTube| = "
                  << std::scientific << std::setprecision(3) << dmax << std::defaultfloat << "\n";
      SCTL_ASSERT_MSG(dmax < (Real)1e-14, "BuildProfiledTube(const) != BuildTube: static path regressed");
    }
    {
      const Long Ntrans = 1;
      const Real xt = (Real)Ntrans * L / (Real)Nelem;
      const Real A = 0.03, k = 2 * const_pi<Real>() / L;
      auto make_prof = [&](Real R_cap) {
        std::function<Real(Real)> r_wall  = [R_cap, A, k](Real x){ return R_cap + A * std::sin(k*x); };
        std::function<Real(Real)> rd_wall = [](Real){ return (Real)0; };
        return glymphbie::MakeCappedRadiusProfile<Real>(R_cap, xt, xt, L, r_wall, rd_wall).r;
      };
      ReportProfiledWatertight<Real>("inner moving", Nelem, cheb, fourier, capOrd, capNaz,
                                     L, r_inner, y0, z0, make_prof(r_inner), tol, comm);
      ReportProfiledWatertight<Real>("outer moving", Nelem, cheb, fourier, capOrd, capNaz,
                                     L, r_outer, y0, z0, make_prof(r_outer), tol, comm);
    }

    // Optional VTK for visual inspection (gitignored vis/ dir convention).
    inner_tube.WriteVTK("vis/capsule_inner_tube", Vector<Real>(), comm);
    inner_caps.WriteVTK("vis/capsule_inner_caps", Vector<Real>(), comm);
    outer_tube.WriteVTK("vis/capsule_outer_tube", Vector<Real>(), comm);
    outer_caps.WriteVTK("vis/capsule_outer_caps", Vector<Real>(), comm);
  }
  Comm::MPI_Finalize();
  return 0;
}
