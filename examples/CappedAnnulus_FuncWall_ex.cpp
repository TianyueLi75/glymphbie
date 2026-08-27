/**
 * CappedAnnulus_FuncWall_ex.cpp
 *
 * Couples the wall-function dynamics (FuncWall) to the NON-PERIODIC (capped)
 * annular Stokes solver (include/CappedAnnulusFlow.hpp). It demonstrates the
 * design from the compatibility plan:
 *
 *   [fixed dome cap] [transition panels] [moving wall] [transition panels] [fixed dome cap]
 *
 * The two hemispherical dome caps stay FIXED (constant radius R_in / R_out); a
 * few CSBQ transition panels at each end blend the fixed cap radius into the
 * FuncWall-defined moving-wall radius r(x,t) using the degree-5 partition of
 * unity (pou_s5). The wall motion enters the BVP as a purely radial slip
 * velocity on the tubes (BuildMovingWallSlipBC), no-slip on the fixed caps.
 *
 * The FuncWall here is a peristaltic breathing wave applied identically to the
 * inner and outer walls: r(x,t) = R_base + A sin(2 pi x / L) sin(omega t). Using
 * the SAME wave shape on both walls (and the symmetric transition bands) makes
 * the net wall flux integrate to ~0, i.e. the motion is volume-conserving, so a
 * fully-closed fixed-cap cavity is compatible (interior-Dirichlet
 * incompressibility) with no compensating cap flux -- we assert this.
 *
 * Usage: ./CappedAnnulus_FuncWall_ex [Nelem fourier capOrder capNaz Nsteps A tol]
 *   Slender CSBQ ElemOrder is fixed at 10; capOrder in {4,8,...,48}.
 * Serial build (no PVFMM/MPI); mirrors the other CappedAnnulus_* standalone
 * examples' include order and config.
 */
#include <sctl.hpp>
#include <quad_junctions/collar_mount_geom.hpp>
#include <quad_junctions/quad_scheme.hpp>
#include <csbq/slender_element.hpp>
#include <csbq/slender_element.cpp>
#include <StokesBIO.hpp>
#include <StokesBIO.cpp>            // header-only use of the project StokesBIO (no GLYMPHBIE link)
#include <CappedAnnulusFlow.hpp>
#include <FuncWall.hpp>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace sctl;

// Peristaltic breathing wall: r(x,t) = R_base + A sin(2 pi x / L) sin(omega t),
// with matching rdot. Uses the (radius, rdot, time, coords) FuncWall signature.
template <class Real> struct PeristalticWall {
  Real R_base, A, omega, L;
  PeristalticWall(Real R_base_, Real A_, Real omega_, Real L_)
      : R_base(R_base_), A(A_), omega(omega_), L(L_) {}
  void operator()(Vector<Real>& radius, Vector<Real>& rdot, Real t, Vector<Real>& coords) const {
    const Real k = 2 * const_pi<Real>() / L;
    for (Long i = 0; i < radius.Dim(); i++) {
      const Real x = coords[i*3+0];
      radius[i] = R_base + A * std::sin(k*x) * std::sin(omega*t);
      rdot[i]   =          A * std::sin(k*x) * omega * std::cos(omega*t);
    }
  }
};

// Adapter: turn a FuncWall into the scalar axial profiles r(x)/rdot(x) that
// CappedAnnulusFlow::Config consumes, frozen at time t. Evaluates the wall
// functor at a single node [x, y0, z0] via FuncWall::evaluate{Inner,Outer}.
template <class Real, class Wall>
void MakeWallProfiles(Wall& wall, Real t, Real y0, Real z0,
                      std::function<Real(Real)>& r_in,  std::function<Real(Real)>& rdot_in,
                      std::function<Real(Real)>& r_out, std::function<Real(Real)>& rdot_out) {
  r_in = [&wall, t, y0, z0](Real x) -> Real {
    Vector<Real> c(3); c[0]=x; c[1]=y0; c[2]=z0; Vector<Real> r, rd;
    wall.evaluateInner(c, t, r, rd); return r[0];
  };
  rdot_in = [&wall, t, y0, z0](Real x) -> Real {
    Vector<Real> c(3); c[0]=x; c[1]=y0; c[2]=z0; Vector<Real> r, rd;
    wall.evaluateInner(c, t, r, rd); return rd[0];
  };
  r_out = [&wall, t, y0, z0](Real x) -> Real {
    Vector<Real> c(3); c[0]=x; c[1]=y0; c[2]=z0; Vector<Real> r, rd;
    wall.evaluateOuter(c, t, r, rd); return r[0];
  };
  rdot_out = [&wall, t, y0, z0](Real x) -> Real {
    Vector<Real> c(3); c[0]=x; c[1]=y0; c[2]=z0; Vector<Real> r, rd;
    wall.evaluateOuter(c, t, r, rd); return rd[0];
  };
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  using Real = double;
  {
    const Comm comm = Comm::World();

    glymphbie::CappedAnnulusFlow<Real>::Config cfg;
    cfg.cheb     = 10;
    cfg.Nelem    = (argc > 1) ? std::atol(argv[1]) : 8;
    cfg.fourier  = (argc > 2) ? std::atoi(argv[2]) : 16;
    cfg.capOrder = (argc > 3) ? std::atoi(argv[3]) : 8;
    cfg.capNaz   = (argc > 4) ? std::atoi(argv[4]) : 16;
    const Long Nsteps = (argc > 5) ? std::atol(argv[5]) : 3;
    const Real A      = (argc > 6) ? std::atof(argv[6]) : 0.02;   // breathing amplitude
    cfg.tol = (argc > 7) ? std::atof(argv[7]) : 1e-8;
    cfg.R_in = 0.3; cfg.R_out = 0.48; cfg.L = 1.0;
    cfg.Ntrans = 1;   // one transition panel per tube end

    const Real dt = 0.05, omega = 2 * const_pi<Real>();  // one breathing period ~ 1 time unit

    // Build a FuncWall whose stored grid matches the tube's Chebyshev x-nodes
    // (only used to hold state / advance time; geometry is sampled via the
    // functor adapter). The grid is the annular centerline node layout.
    Vector<Real> grid; for (Long k = 0; k < cfg.Nelem; k++) {
      const Vector<Real>& cn = SlenderElemList<Real>::CenterlineNodes(cfg.cheb);
      for (Integer j = 0; j < cfg.cheb; j++) {
        const Real x = cfg.L * (k + cn[j]) / cfg.Nelem;
        grid.PushBack(x); grid.PushBack(cfg.y0); grid.PushBack(cfg.z0);
      }
    }
    const Long Ng = grid.Dim()/3;
    Vector<Real> r_in0(Ng), r_out0(Ng); r_in0 = cfg.R_in; r_out0 = cfg.R_out;

    PeristalticWall<Real> inner_fn(cfg.R_in,  A, omega, cfg.L);
    PeristalticWall<Real> outer_fn(cfg.R_out, A, omega, cfg.L);
    FuncWall<Real, PeristalticWall<Real>, PeristalticWall<Real>>
        wall(dt, grid, grid, r_out0, r_in0, inner_fn, outer_fn);

    if (!comm.Rank())
      std::cout << "Capped-annulus + FuncWall (Nelem=" << cfg.Nelem << " fourier=" << cfg.fourier
                << " capOrder=" << cfg.capOrder << " Ntrans=" << cfg.Ntrans
                << " A=" << A << " Nsteps=" << Nsteps << ")\n"
                << "  R_in=" << cfg.R_in << " R_out=" << cfg.R_out << " L=" << cfg.L << "\n";

    for (Long step = 0; step < Nsteps; step++) {
      const Real t = wall.getTime();

      // Freeze the wall-function profiles at the current time and (re)build the
      // capped geometry: fixed caps + blended transition + moving wall.
      MakeWallProfiles<Real>(wall, t, cfg.y0, cfg.z0,
                             cfg.r_in_of_x, cfg.rdot_in_of_x, cfg.r_out_of_x, cfg.rdot_out_of_x);
      glymphbie::CappedAnnulusFlow<Real> flow(cfg, comm);

      // Moving-wall slip BC; check volume conservation (net tube flux ~ 0).
      Real net_flux = 0;
      Vector<Real> bc = flow.BuildMovingWallSlipBC(&net_flux);
      if (!comm.Rank())
        std::cout << std::scientific << std::setprecision(3)
                  << "  step " << step << " t=" << std::fixed << std::setprecision(3) << t
                  << std::scientific << "  net wall flux = " << net_flux
                  << "  (want ~0: volume preservation)\n" << std::defaultfloat;
      SCTL_ASSERT_MSG(std::fabs(net_flux) < 1e-6,
                      "net wall flux not ~0: wall motion is not volume-conserving; "
                      "add a compensating cap flux (BuildInflowOutflowBC)");

      // Solve and sample the interior axial velocity at the mid-section.
      Vector<Real> sigma = flow.Solve(bc, 400, "capped-funcwall");

      const Long Nr = 5, Nth = 8;
      Vector<Real> Xtrg;
      for (Long i = 0; i < Nr; i++) {
        const Real r = cfg.R_in + (cfg.R_out - cfg.R_in) * (Real)(i+1) / (Real)(Nr+1);
        for (Long j = 0; j < Nth; j++) {
          const Real th = 2*const_pi<Real>()*(Real)j/(Real)Nth;
          Xtrg.PushBack(cfg.L/2);
          Xtrg.PushBack(cfg.y0 + r*std::cos(th));
          Xtrg.PushBack(cfg.z0 + r*std::sin(th));
        }
      }
      Vector<Real> Utrg; flow.EvaluateVelocity(sigma, Xtrg, Utrg);
      Real umax = 0; for (Long i = 0; i < Utrg.Dim()/3; i++) {
        const Real u2 = Utrg[i*3]*Utrg[i*3] + Utrg[i*3+1]*Utrg[i*3+1] + Utrg[i*3+2]*Utrg[i*3+2];
        umax = std::max(umax, std::sqrt(u2));
      }
      if (!comm.Rank())
        std::cout << std::scientific << std::setprecision(4)
                  << "    mid-section max |u| = " << umax << std::defaultfloat << "\n";

      wall.update();   // advance wall state + time to the next step
    }
  }
  Comm::MPI_Finalize();
  return 0;
}
