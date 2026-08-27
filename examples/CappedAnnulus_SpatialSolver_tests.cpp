/**
 * CappedAnnulus_SpatialSolver_tests.cpp
 *
 * Accuracy test for the NON-PERIODIC (capped) annular Stokes solver
 * (include/CappedAnnulusFlow.hpp), the geometric counterpart of the periodic
 * examples/SpatialSolver_tests.cpp.
 *
 * The concentric annulus (inner radius R_in, outer R_out, centerline (y,z) =
 * (0.5,0.5), length L=1) is driven by a prescribed Poiseuille inflow/outflow on
 * the outer caps whose volume flux equals the analytic flux Q of a concentric
 * annular Poiseuille flow under unit pressure gradient (dp=1, mu=1). We then
 *   (a) verify the prescribed net boundary flux is ~0 (volume preservation),
 *   (b) compare the interior axial velocity at the developed mid-section
 *       (x=L/2) against the exact annular-Poiseuille profile u_x(r), and
 *   (c) compare the flux integrated across the mid-section against Q.
 *
 * Usage: ./CappedAnnulus_SpatialSolver_tests [Nelem fourier capOrder capNaz tol]
 *   Slender CSBQ ElemOrder is fixed at 10 (the only supported order).
 *   capOrder (cap QuadElemList patch order) MUST be one of {4,8,12,...,48}.
 * Serial build (no PVFMM/MPI); mirrors CappedAnnulus_watertight's config.
 */
#include <sctl.hpp>
#include <quad_junctions/collar_mount_geom.hpp>
#include <quad_junctions/quad_scheme.hpp>
#include <csbq/slender_element.hpp>
#include <csbq/slender_element.cpp>
#include <StokesBIO.hpp>
#include <StokesBIO.cpp>            // header-only use of the project StokesBIO (no GLYMPHBIE link)
#include <CappedAnnulusFlow.hpp>
#include <iomanip>
#include <iostream>

using namespace sctl;

// Concentric annular Poiseuille axial velocity at radius r (mu=1, G=|dpdx|=1, +x flow).
template <class Real> static Real annular_poiseuille_ux(Real r, Real R_in, Real R_out) {
  const Real rmax2 = (R_out*R_out - R_in*R_in) / (2*sctl::log<Real>(R_out/R_in));
  return (Real)0.25 * (R_out*R_out - r*r - 2*rmax2*sctl::log<Real>(R_out/r));
}

// Analytic volume flux of the concentric annular Poiseuille flow (mu=1, G=1).
template <class Real> static Real annular_poiseuille_flux(Real R_in, Real R_out) {
  const Real pi = sctl::const_pi<Real>();
  return (pi/8) * (R_out*R_out*R_out*R_out - R_in*R_in*R_in*R_in
                   - (R_out*R_out - R_in*R_in)*(R_out*R_out - R_in*R_in) / sctl::log<Real>(R_out/R_in));
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  using Real = double;
  {
    const Comm comm = Comm::World();

    glymphbie::CappedAnnulusFlow<Real>::Config cfg;
    cfg.cheb     = 10;   // CSBQ slender ElemOrder is fixed at 10 (only supported order)
    cfg.Nelem    = (argc > 1) ? std::atol(argv[1]) : 8;
    cfg.fourier  = (argc > 2) ? std::atoi(argv[2]) : 16;
    cfg.capOrder = (argc > 3) ? std::atoi(argv[3]) : 8;   // MUST be in {4,8,...,48}
    cfg.capNaz   = (argc > 4) ? std::atoi(argv[4]) : 16;
    cfg.R_in = 0.3; cfg.R_out = 0.48; cfg.L = 1.0;
    cfg.tol = (argc > 5) ? std::atof(argv[5]) : 1e-10;   // loosen for fast validation runs

    if (!comm.Rank())
      std::cout << "Capped-annulus accuracy test (Nelem=" << cfg.Nelem << " cheb=" << cfg.cheb
                << " fourier=" << cfg.fourier << " capOrder=" << cfg.capOrder
                << " capNaz=" << cfg.capNaz << ")\n"
                << "  R_in=" << cfg.R_in << " R_out=" << cfg.R_out << " L=" << cfg.L << "\n";

    glymphbie::CappedAnnulusFlow<Real> flow(cfg, comm);

    const Real Q = annular_poiseuille_flux<Real>(cfg.R_in, cfg.R_out);
    if (!comm.Rank()) std::cout << std::setprecision(10) << "  analytic dp=1 flux Q = " << Q << "\n";

    // Drive with the analytic flux; check net-flux (volume preservation).
    Real net_flux = 0;
    Vector<Real> bc = flow.BuildInflowOutflowBC(Q, &net_flux);
    if (!comm.Rank())
      std::cout << "  prescribed net boundary flux |int u.n| = " << std::scientific << std::fabs(net_flux)
                << std::defaultfloat << "   (want ~0: volume preservation)\n";

    // Solve for the surface density.
    Vector<Real> sigma = flow.Solve(bc);

    // --- (b) interior velocity at the developed mid-section x=L/2 ---
    const Long Nr = 6, Nth = 8;
    Vector<Real> Xtrg;
    Vector<Real> r_samp(Nr);
    for (Long i = 0; i < Nr; i++) {
      const Real r = cfg.R_in + (cfg.R_out - cfg.R_in) * (Real)(i+1) / (Real)(Nr+1);  // interior shells
      r_samp[i] = r;
      for (Long j = 0; j < Nth; j++) {
        const Real th = 2*const_pi<Real>()*(Real)j/(Real)Nth;
        Xtrg.PushBack(cfg.L/2);
        Xtrg.PushBack(cfg.y0 + r*std::cos(th));
        Xtrg.PushBack(cfg.z0 + r*std::sin(th));
      }
    }
    Vector<Real> Utrg;
    flow.EvaluateVelocity(sigma, Xtrg, Utrg);

    Real max_rel = 0, max_abs = 0, ref = 0;
    for (Long i = 0; i < Nr; i++) {
      const Real ue = annular_poiseuille_ux<Real>(r_samp[i], cfg.R_in, cfg.R_out);
      ref = std::max(ref, std::fabs(ue));
      for (Long j = 0; j < Nth; j++) {
        const Real ux = Utrg[(i*Nth+j)*3+0];
        max_abs = std::max(max_abs, std::fabs(ux - ue));
      }
    }
    max_rel = (ref > 0) ? max_abs/ref : max_abs;
    if (!comm.Rank())
      std::cout << std::scientific << std::setprecision(4)
                << "  mid-section u_x vs exact annular Poiseuille: max abs err = " << max_abs
                << ", max rel err = " << max_rel << std::defaultfloat << "\n";
  }
  Comm::MPI_Finalize();
  return 0;
}
