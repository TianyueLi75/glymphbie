/**
 * CappedAnnulus_TracerParticles_verify.cpp
 *
 * Verification that the tracer-particle velocity sampling (TracerParticles +
 * CappedAnnulusFlow::EvaluateVelocity) is correct, on the simplest case: a
 * STATIC, constant-radius concentric annulus (drdt = 0, a "straight channel")
 * driven to the exact annular Poiseuille flow.
 *
 * The concentric annulus (R_in, R_out, centerline (0.5,0.5), L=1) is driven by
 * the prescribed outer-cap inflow/outflow whose volume flux equals the analytic
 * annular-Poiseuille flux Q (dp=1, mu=1) -- exactly as in
 * CappedAnnulus_SpatialSolver_tests.cpp. We then seed a handful of tracers at
 * known radii on the developed mid-section (x=L/2), evaluate the flow velocity
 * AT the tracer positions through the tracer's own VelFn wrapper, and assert the
 * axial component matches the exact profile u_x(r). Finally we take one explicit
 * TracerParticles::Step and assert the displacement equals dt * u (the forward
 * Euler integrator path).
 *
 * Usage: ./CappedAnnulus_TracerParticles_verify [Nelem fourier capOrder capNaz tol]
 *   Slender CSBQ ElemOrder is fixed at 10; capOrder in {4,8,...,48}.
 * Serial build (no PVFMM/MPI).
 */
#include <sctl.hpp>
#include <quad_junctions/collar_mount_geom.hpp>
#include <quad_junctions/quad_scheme.hpp>
#include <csbq/slender_element.hpp>
#include <csbq/slender_element.cpp>
#include <StokesBIO.hpp>
#include <StokesBIO.cpp>            // header-only use of the project StokesBIO (no GLYMPHBIE link)
#include <CappedAnnulusFlow.hpp>
#include <Particle.hpp>
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
  int rc = 0;
  {
    const Comm comm = Comm::World();

    glymphbie::CappedAnnulusFlow<Real>::Config cfg;
    cfg.cheb     = 10;   // CSBQ slender ElemOrder is fixed at 10 (only supported order)
    cfg.Nelem    = (argc > 1) ? std::atol(argv[1]) : 6;
    cfg.fourier  = (argc > 2) ? std::atoi(argv[2]) : 16;
    cfg.capOrder = (argc > 3) ? std::atoi(argv[3]) : 8;   // MUST be in {4,8,...,48}
    cfg.capNaz   = (argc > 4) ? std::atoi(argv[4]) : 16;
    cfg.R_in = 0.3; cfg.R_out = 0.48; cfg.L = 1.0;
    cfg.tol = (argc > 5) ? std::atof(argv[5]) : 1e-8;
    // Static geometry: leave all Config::*_of_x empty => constant radius, drdt=0.

    if (!comm.Rank())
      std::cout << "Tracer-particle velocity verification (static drdt=0 annular Poiseuille)\n"
                << "  Nelem=" << cfg.Nelem << " cheb=" << cfg.cheb << " fourier=" << cfg.fourier
                << " capOrder=" << cfg.capOrder << " capNaz=" << cfg.capNaz << "\n"
                << "  R_in=" << cfg.R_in << " R_out=" << cfg.R_out << " L=" << cfg.L << "\n";

    glymphbie::CappedAnnulusFlow<Real> flow(cfg, comm);

    // Drive with the analytic annular-Poiseuille flux; check net-flux (volume preservation).
    const Real Q = annular_poiseuille_flux<Real>(cfg.R_in, cfg.R_out);
    Real net_flux = 0;
    Vector<Real> bc = flow.BuildInflowOutflowBC(Q, &net_flux);
    if (!comm.Rank())
      std::cout << std::setprecision(10) << "  analytic dp=1 flux Q = " << Q
                << ",  prescribed |net boundary flux| = " << std::scientific << std::fabs(net_flux)
                << std::defaultfloat << "  (want ~0)\n";

    Vector<Real> sigma = flow.Solve(bc);

    // Seed tracers at known radii on the developed mid-section x=L/2.
    const Long Np = 6;
    Vector<Real> r_samp(Np), X0;
    for (Long i = 0; i < Np; i++) {
      const Real r = cfg.R_in + (cfg.R_out - cfg.R_in) * (Real)(i+1) / (Real)(Np+1);
      r_samp[i] = r;
      X0.PushBack(cfg.L/2);           // x
      X0.PushBack(cfg.y0 + r);        // y (theta = 0)
      X0.PushBack(cfg.z0);            // z
    }

    typename TracerParticles<Real>::Geometry geom;
    geom.x0 = 0; geom.x1 = cfg.L; geom.y0 = cfg.y0; geom.z0 = cfg.z0;
    geom.R_in = cfg.R_in; geom.R_out = cfg.R_out;
    TracerParticles<Real> tracers(X0, geom);

    // The tracer's velocity evaluator: wrap CappedAnnulusFlow::EvaluateVelocity.
    auto vel = [&](const Vector<Real>& X, Vector<Real>& U) { flow.EvaluateVelocity(sigma, X, U); };

    // (a) velocity sampled at tracer positions vs exact annular Poiseuille u_x.
    Vector<Real> U; vel(tracers.Positions(), U);
    Real max_abs = 0, ref = 0;
    for (Long i = 0; i < Np; i++) {
      const Real ue = annular_poiseuille_ux<Real>(r_samp[i], cfg.R_in, cfg.R_out);
      ref = std::max(ref, std::fabs(ue));
      max_abs = std::max(max_abs, std::fabs(U[i*3+0] - ue));
    }
    const Real max_rel = (ref > 0) ? max_abs/ref : max_abs;
    if (!comm.Rank())
      std::cout << std::scientific << std::setprecision(4)
                << "  tracer-sampled u_x vs exact: max abs err = " << max_abs
                << ", max rel err = " << max_rel << std::defaultfloat << "\n";

    // (b) one explicit Step advances positions by dt*u (forward Euler path).
    const Real dt = 0.01;
    Vector<Real> Xbefore = tracers.Positions();
    tracers.Step(dt, vel);
    Real max_step_err = 0;
    const Vector<Real>& Xafter = tracers.Positions();
    for (Long i = 0; i < Xafter.Dim(); i++)
      max_step_err = std::max(max_step_err, std::fabs(Xafter[i] - (Xbefore[i] + dt*U[i])));
    if (!comm.Rank())
      std::cout << std::scientific << std::setprecision(4)
                << "  Step displacement vs dt*u: max err = " << max_step_err << std::defaultfloat << "\n";

    // Pass/fail: the sampled velocity should match to solver accuracy; the Euler
    // step must match dt*u to round-off (no re-seed at the mid-section).
    const Real vel_tol = 1e-2;   // BIE @ tol~1e-8 with a small Nelem: a loose but meaningful bound
    const bool ok = (max_rel < vel_tol) && (max_step_err < 1e-12);
    if (!comm.Rank())
      std::cout << (ok ? "  PASS" : "  FAIL")
                << " (vel_rel<" << vel_tol << " && step_err<1e-12)\n";
    rc = ok ? 0 : 1;
  }
  Comm::MPI_Finalize();
  return rc;
}
