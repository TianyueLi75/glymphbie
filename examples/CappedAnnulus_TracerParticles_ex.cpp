/**
 * CappedAnnulus_TracerParticles_ex.cpp
 *
 * End-to-end wall-function-driven capped-annulus flow with Lagrangian tracer
 * particles and a flow-field VTK time series. Combines:
 *   - the wall-function coupling of CappedAnnulus_FuncWall_ex.cpp (a peristaltic
 *     FuncWall frozen into CappedAnnulusFlow radius/rdot profiles each step),
 *   - a SUPERIMPOSED inflow/outflow (pin-pout) flux Q on the outer caps, so the
 *     flow is driven by BOTH the moving-wall slip AND a compatible pressure-like
 *     throughput. The two BCs have disjoint support (slip on the tubes, flux on
 *     the outer caps), so their AoS BC vectors are summed exactly; the combined
 *     net flux |net_wall + net_cap| must be ~0 (interior-Dirichlet
 *     incompressibility / volume-preservation compatibility), which we assert,
 *   - TracerParticles (Particle.hpp) advected by the solved velocity field, and
 *   - a small interior flow-field grid (FlowVis.hpp) written each step.
 *
 * Output (into vis/): exactly TWO files, each a single VTKHDF (.vtkhdf) time
 * series that ParaView opens natively (no per-step .vtu/.pvtu/.pvd clutter):
 *   - vis/walls.vtkhdf   : the changing wall geometry every step (inner/outer
 *     tubes + caps merged into one grid, colored by the Dirichlet BC and an
 *     ObjectId array: 0=inner tube, 1=outer tube, 2=inner caps, 3=outer caps),
 *   - vis/tracers.vtkhdf : the tracer point cloud every step, colored by the
 *     local velocity (this IS the time-series/trajectory data -- use ParaView's
 *     "Temporal Particles To Pathlines" to draw trajectories from it).
 * The per-class .vtu/.pvd writers are retained on the classes but no longer
 * called here, so a run produces one downloadable file per thing regardless of
 * the number of time steps.
 *
 * Usage: ./CappedAnnulus_TracerParticles_ex [Nelem fourier capOrder capNaz Nsteps A Q tol]
 *   Slender CSBQ ElemOrder is fixed at 10; capOrder in {4,8,...,48}. Defaults are
 *   deliberately small (verification-scale) grids/steps.
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
#include <FuncWall.hpp>
#include <Particle.hpp>
#include <FlowVis.hpp>
#include <VtkHdfWriter.hpp>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

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
    cfg.Nelem    = (argc > 1) ? std::atol(argv[1]) : 6;
    cfg.fourier  = (argc > 2) ? std::atoi(argv[2]) : 16;
    cfg.capOrder = (argc > 3) ? std::atoi(argv[3]) : 8;
    cfg.capNaz   = (argc > 4) ? std::atoi(argv[4]) : 16;
    const Long Nsteps = (argc > 5) ? std::atol(argv[5]) : 4;
    const Real A      = (argc > 6) ? std::atof(argv[6]) : 0.02;   // breathing amplitude
    const Real Q      = (argc > 7) ? std::atof(argv[7]) : 0.0;    // pin-pout throughput flux
    cfg.tol = (argc > 8) ? std::atof(argv[8]) : 1e-8;
    cfg.R_in = 0.3; cfg.R_out = 0.48; cfg.L = 1.0;
    cfg.Ntrans = 1;   // one transition panel per tube end

    const Real dt = 0.05, omega = 2 * const_pi<Real>();  // one breathing period ~ 1 time unit

    // Tracer seeding resolution: tracers are seeded on a uniform (axial, radial,
    // azimuthal) lattice filling the annular gap. THESE are the variables that
    // set how many tracers there are and how densely they fill the domain --
    // raise them for more, more uniformly distributed tracers.
    const Long Nseed_x = 5, Nseed_r = 4, Nseed_th = 8;
    const Long Np = Nseed_x * Nseed_r * Nseed_th;   // total tracer count

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

    // Tracer geometry (fixed cap radii; adequate for the small verification run).
    typename TracerParticles<Real>::Geometry geom;
    geom.x0 = 0; geom.x1 = cfg.L; geom.y0 = cfg.y0; geom.z0 = cfg.z0;
    geom.R_in = cfg.R_in; geom.R_out = cfg.R_out;

    // Seed tracers on a uniform (x, r, theta) lattice filling the annular gap:
    // axial cell centers in (0,L), radial cell centers in (R_in,R_out), azimuth
    // evenly around the annulus. Cell-centered offsets keep seeds off the walls.
    Vector<Real> Xp;
    for (Long ix = 0; ix < Nseed_x; ix++)
      for (Long ir = 0; ir < Nseed_r; ir++)
        for (Long it = 0; it < Nseed_th; it++) {
          const Real x  = cfg.L * (ix + (Real)0.5) / Nseed_x;
          const Real r  = cfg.R_in + (cfg.R_out - cfg.R_in) * (ir + (Real)0.5) / Nseed_r;
          const Real th = 2*const_pi<Real>() * it / Nseed_th;
          Xp.PushBack(x);
          Xp.PushBack(cfg.y0 + r*std::cos(th));
          Xp.PushBack(cfg.z0 + r*std::sin(th));
        }
    TracerParticles<Real> tracers(Xp, geom);

    if (!comm.Rank())
      std::cout << "Capped-annulus + FuncWall + tracers (Nelem=" << cfg.Nelem << " fourier=" << cfg.fourier
                << " capOrder=" << cfg.capOrder << " Ntrans=" << cfg.Ntrans
                << " A=" << A << " Q=" << Q << " Nsteps=" << Nsteps << ")\n"
                << "  R_in=" << cfg.R_in << " R_out=" << cfg.R_out << " L=" << cfg.L
                << ", " << Np << " tracers\n";

    // One VTKHDF time series per "thing": the moving wall geometry and the
    // tracer cloud. Each accumulates all steps in memory and is written as a
    // single .vtkhdf file at the end (see VtkHdfWriter.hpp).
    glymphbie::VtkHdfTimeSeries walls_ts, tracers_ts;

    // Append the current wall surfaces (colored by the split Dirichlet BC) as one
    // merged grid step, and the current tracer cloud (colored by velocity Up) as
    // one point-cloud step. ObjectIds tag the four wall blocks so ParaView can
    // split/color them; the caps ride along every step (self-contained snapshots).
    auto push_snapshots = [&](Real t, glymphbie::CappedAnnulusFlow<Real>& flow,
                              const Vector<Real>& bc, const Vector<Real>& Up) {
      Vector<Real> F_it, F_ot, F_ic, F_oc;
      flow.SplitSurfaceField(bc, F_it, F_ot, F_ic, F_oc);
      sctl::VTUData v_it, v_ot, v_ic, v_oc;
      flow.GetSurfaceVTUData(v_it, v_ot, v_ic, v_oc, F_it, F_ot, F_ic, F_oc);
      glymphbie::VtkHdfStep wstep;
      wstep.AppendVTU(v_it, "BC", 0);
      wstep.AppendVTU(v_ot, "BC", 1);
      wstep.AppendVTU(v_ic, "BC", 2);
      wstep.AppendVTU(v_oc, "BC", 3);
      walls_ts.AddStep((double)t, std::move(wstep));

      sctl::VTUData v_tr;
      tracers.GetVTUData(v_tr, Up);
      glymphbie::VtkHdfStep tstep;
      tstep.AppendVTU(v_tr, "Velocity", 0);
      tracers_ts.AddStep((double)t, std::move(tstep));
    };

    for (Long step = 0; step < Nsteps; step++) {
      const Real t = wall.getTime();

      // Freeze the wall-function profiles at the current time and (re)build the
      // capped geometry: fixed caps + blended transition + moving wall.
      MakeWallProfiles<Real>(wall, t, cfg.y0, cfg.z0,
                             cfg.r_in_of_x, cfg.rdot_in_of_x, cfg.r_out_of_x, cfg.rdot_out_of_x);
      glymphbie::CappedAnnulusFlow<Real> flow(cfg, comm);

      // Combined Dirichlet BC: moving-wall slip (tubes) + inflow/outflow flux Q
      // (outer caps). Disjoint support => sum is exact. Check compatibility.
      Real net_wall = 0, net_cap = 0;
      Vector<Real> bc = flow.BuildMovingWallSlipBC(&net_wall);
      if (Q != (Real)0) {
        Vector<Real> bc_cap = flow.BuildInflowOutflowBC(Q, &net_cap);
        for (Long i = 0; i < bc.Dim(); i++) bc[i] += bc_cap[i];
      }
      const Real net_total = net_wall + net_cap;
      if (!comm.Rank())
        std::cout << std::scientific << std::setprecision(3)
                  << "  step " << step << " t=" << std::fixed << std::setprecision(3) << t
                  << std::scientific << "  net flux: wall=" << net_wall << " cap=" << net_cap
                  << " total=" << net_total << "  (want total ~0)\n" << std::defaultfloat;
      SCTL_ASSERT_MSG(std::fabs(net_total) < 1e-6,
                      "combined net flux not ~0: wall motion + pin-pout flux are not "
                      "volume-compatible for a fixed-cap interior-Dirichlet cavity");

      Vector<Real> sigma = flow.Solve(bc, 400, "capped-tracers");

      // The tracer / field velocity evaluator wraps EvaluateVelocity.
      auto vel = [&](const Vector<Real>& X, Vector<Real>& U) { flow.EvaluateVelocity(sigma, X, U); };

      // Accumulate this step into the two VTKHDF time series: the wall geometry
      // colored by the split Dirichlet BC (moving-wall slip on the tubes +
      // inflow/outflow on the outer caps), and the tracer cloud colored by the
      // current velocity.
      Vector<Real> Up; vel(tracers.Positions(), Up);
      push_snapshots(t, flow, bc, Up);

      // Record history, then advect.
      tracers.RecordHistory(t);
      tracers.Step(dt, vel);   // forward Euler with the current (frozen) velocity field

      wall.update();           // advance wall state + time to the next step
    }

    // Final snapshot at the final wall state. wall.update() already advanced the
    // wall Nsteps times in the loop, so getTime() is the endpoint; rebuild the
    // flow there (the loop's `flow`/`vel` are out of scope) for one last geometry
    // + tracer + velocity snapshot, matched to the final tracer positions.
    {
      const Real t = wall.getTime();
      MakeWallProfiles<Real>(wall, t, cfg.y0, cfg.z0,
                             cfg.r_in_of_x, cfg.rdot_in_of_x, cfg.r_out_of_x, cfg.rdot_out_of_x);
      glymphbie::CappedAnnulusFlow<Real> flow(cfg, comm);

      Real net_wall = 0, net_cap = 0;
      Vector<Real> bc = flow.BuildMovingWallSlipBC(&net_wall);
      if (Q != (Real)0) {
        Vector<Real> bc_cap = flow.BuildInflowOutflowBC(Q, &net_cap);
        for (Long i = 0; i < bc.Dim(); i++) bc[i] += bc_cap[i];
      }
      Vector<Real> sigma = flow.Solve(bc, 400, "capped-tracers-final");
      auto vel = [&](const Vector<Real>& X, Vector<Real>& U) { flow.EvaluateVelocity(sigma, X, U); };

      // Final recorded position, then the final wall + tracer snapshot.
      tracers.RecordHistory(t);
      Vector<Real> Up; vel(tracers.Positions(), Up);
      push_snapshots(t, flow, bc, Up);
    }

    // Write the two accumulated time series -- ONE file each, however many steps.
    if (!comm.Rank()) {
      std::filesystem::create_directories("vis");
      bool ok_w = walls_ts.Write("vis/walls");
      bool ok_t = tracers_ts.Write("vis/tracers");
      if (ok_w && ok_t) {
        std::cout << "  wrote vis/walls.vtkhdf (" << walls_ts.NumSteps() << " steps) and "
                     "vis/tracers.vtkhdf (" << tracers_ts.NumSteps() << " steps)\n";
      } else {
        std::cerr << "  ERROR: failed to write "
                  << (!ok_w ? "vis/walls.vtkhdf " : "")
                  << (!ok_t ? "vis/tracers.vtkhdf" : "") << "\n";
        return 1;
      }
    }
  }
  Comm::MPI_Finalize();
  return 0;
}
