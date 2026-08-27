#ifndef _CAPPEDANNULUSFLOW_HPP_
#define _CAPPEDANNULUSFLOW_HPP_

/**
 * CappedAnnulusFlow.hpp
 *
 * Non-periodic ("capped") Stokes flow in the annular perivascular space. Where
 * the periodic framing (StokesBIO::SetPeriodicity + a background pressure
 * gradient) closes the channel by periodicity in x, this closes it
 * geometrically: the inner (arteriole) and outer (glial) walls are slender
 * tubes, each sealed at both ends by a hemispherical butterfly dome cap
 * (quad-junctions add_tip_cap_butterfly), giving two nested watertight capsules
 * whose gap is the fluid cavity. See include/CappedAnnulusGeom.hpp for the
 * geometry builders and examples/CappedAnnulus_watertight.cpp for the
 * closure/watertightness checks this is built on.
 *
 * The fluid domain is doubly connected: interior to the OUTER capsule and
 * exterior to the INNER capsule. We solve the interior Stokes velocity
 * (Dirichlet) BVP with a combined-field representation
 *     u(x) = SL_scal * S[sigma](x) + DL_scal * D[sigma](x),
 * assembled over all four element lists with a single generalized StokesBIO
 * (include/StokesBIO.hpp, which now accepts any sctl::ElementListBase). The DL
 * self-jump is taken from the known outward orientation, per block:
 *   - OUTER wall (tube + caps): CSBQ normal points OUT of the fluid  -> +jump
 *   - INNER wall (tube + caps): CSBQ normal points INTO the fluid    -> -jump
 * (the inner capsule behaves like an obstacle inside the outer container, cf.
 * quad-junctions solve_dirichlet_bvp's obstacle handling). interior => jump =
 * -0.5*DL_scal; we use SL_scal=-1, DL_scal=+1 so SL_scal*DL_scal<0 and the
 * interior CFIE null space is already avoided (no SL sign flip needed).
 *
 * Flow is driven by prescribed inflow/outflow on the OUTER caps only (no-slip on
 * both tubes and the inner caps -- the solid arteriole ends do not pass fluid):
 * a Poiseuille-like axial profile on the low-x outer dome (inflow) and high-x
 * outer dome (outflow), each amplitude-normalized so its cross-cap volume flux
 * equals the prescribed Q. Since inflow flux = outflow flux, the net boundary
 * flux integral(u.n) dA = 0 -- the incompressibility compatibility condition,
 * i.e. VOLUME PRESERVATION -- which we assert.
 *
 * Recipe ported from quad-junctions solve_dirichlet_bvp (hybrid_bie_tests.hpp)
 * and the cap-BC block of src/ybifurc-flow-bie.cpp, generalized from two element
 * lists to the four of the nested-capsule annulus.
 *
 * Header-only; template on Real. Uses the project's generalized StokesBIO, so a
 * driver may build it standalone by including StokesBIO.cpp (as the examples do).
 */

#include <sctl.hpp>
#include <StokesBIO.hpp>
#include "CappedAnnulusGeom.hpp"
#include <string>
#include <cmath>
#include <functional>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

namespace glymphbie {

template <class Real>
class CappedAnnulusFlow {
 public:
  // Discretization + geometry parameters.
  struct Config {
    sctl::Long    Nelem   = 8;     // slender panels per tube
    sctl::Integer cheb    = 10;    // Chebyshev order per slender panel (any order; CSBQ)
    sctl::Integer fourier = 16;    // Fourier order around each tube
    sctl::Integer capOrder = 8;    // cap (QuadElemList) patch order -- MUST be in {4,8,...,48}
    sctl::Integer capNaz  = 16;    // azimuthal panels around each cap equator
    Real L      = 1.0;             // axial length
    Real y0     = 0.5, z0 = 0.5;   // annulus centerline (y,z); matches Annular convention
    Real R_in   = 0.3;             // inner (arteriole) wall radius
    Real R_out  = 0.48;            // outer (glial) wall radius
    Real tol    = 1e-10;           // operator / quadrature accuracy (drives Duffy self/near cost)
    // Cap singular-quadrature knobs for QuadElemList::SetQuadScheme (Duffy default
    // ignores cov_q/Nbeta -- its order comes from tol -- but honours max_depth).
    sctl::Integer cov_q     = 6;   // one of {6,10}
    sctl::Integer Nbeta     = 48;  // one of {48,100,200,300,400,512}
    sctl::Integer max_depth = 4;   // one of {4,8,12,30}

    // --- Moving-wall coupling (optional) --------------------------------------
    // When the *_of_x radius profiles are set, the inner/outer TUBES are built
    // with an axial radius profile that blends the fixed cap radius R_in/R_out
    // into a wall-function-defined radius over Ntrans transition panels at each
    // end (see MakeCappedRadiusProfile / BuildProfiledTube). The dome caps stay
    // fixed at R_in/R_out. Leave the profiles empty (default) for the original
    // constant-radius, statically-capped geometry.
    sctl::Long Ntrans = 1;                    // transition panels per tube end
    std::function<Real(Real)> r_in_of_x;      // inner wall radius r(x); empty => constant R_in
    std::function<Real(Real)> r_out_of_x;     // outer wall radius r(x); empty => constant R_out
    std::function<Real(Real)> rdot_in_of_x;   // inner wall rdot(x); empty => 0 (static)
    std::function<Real(Real)> rdot_out_of_x;  // outer wall rdot(x); empty => 0 (static)
  };

  CappedAnnulusFlow(const Config& cfg, const sctl::Comm& comm)
      : cfg_(cfg), comm_(comm) { Setup(); }

  // (Re)build the four element lists and the combined node/jump-sign data.
  void Setup() {
    // CSBQ slender elements are only used at Chebyshev ElemOrder 10 in this
    // project (the CSBQ singular/near quadrature tables and convergence are
    // validated for order 10; other orders are unsupported here).
    SCTL_ASSERT_MSG(cfg_.cheb == 10, "CappedAnnulusFlow: CSBQ slender ElemOrder (cheb) must be 10");
    // Cap patch order is constrained by SCTL's templated QuadElemList near/self
    // quadrature schemes to {4,8,...,48} (independent of the slender order).
    SCTL_ASSERT_MSG(cfg_.capOrder >= 4 && cfg_.capOrder <= 48 && cfg_.capOrder % 4 == 0,
                    "CappedAnnulusFlow: capOrder must be one of {4,8,12,...,48}");
    // Tubes: either constant-radius (legacy) or a fixed-cap -> moving-wall
    // blended profile. The blended rdot(x) profiles are cached for the
    // moving-wall slip BC (BuildMovingWallSlipBC).
    rdot_in_prof_ = {}; rdot_out_prof_ = {};
    if (cfg_.r_in_of_x || cfg_.r_out_of_x) {
      const Real xt = (Real)cfg_.Ntrans * cfg_.L / (Real)cfg_.Nelem;   // transition band length per end
      auto zero = [](Real) -> Real { return (Real)0; };
      auto r_in_wall    = cfg_.r_in_of_x    ? cfg_.r_in_of_x    : std::function<Real(Real)>([R=cfg_.R_in ](Real){ return R; });
      auto r_out_wall   = cfg_.r_out_of_x   ? cfg_.r_out_of_x   : std::function<Real(Real)>([R=cfg_.R_out](Real){ return R; });
      auto rdot_in_wall  = cfg_.rdot_in_of_x  ? cfg_.rdot_in_of_x  : std::function<Real(Real)>(zero);
      auto rdot_out_wall = cfg_.rdot_out_of_x ? cfg_.rdot_out_of_x : std::function<Real(Real)>(zero);
      auto prof_in  = MakeCappedRadiusProfile<Real>(cfg_.R_in,  xt, xt, cfg_.L, r_in_wall,  rdot_in_wall);
      auto prof_out = MakeCappedRadiusProfile<Real>(cfg_.R_out, xt, xt, cfg_.L, r_out_wall, rdot_out_wall);
      inner_tube_ = BuildProfiledTube<Real>(cfg_.Nelem, cfg_.cheb, cfg_.fourier, cfg_.L, cfg_.y0, cfg_.z0, prof_in.r);
      outer_tube_ = BuildProfiledTube<Real>(cfg_.Nelem, cfg_.cheb, cfg_.fourier, cfg_.L, cfg_.y0, cfg_.z0, prof_out.r);
      rdot_in_prof_  = prof_in.rdot;
      rdot_out_prof_ = prof_out.rdot;
    } else {
      inner_tube_ = BuildTube<Real>(cfg_.Nelem, cfg_.cheb, cfg_.fourier, cfg_.L, cfg_.R_in, cfg_.y0, cfg_.z0);
      outer_tube_ = BuildTube<Real>(cfg_.Nelem, cfg_.cheb, cfg_.fourier, cfg_.L, cfg_.R_out, cfg_.y0, cfg_.z0);
    }
    inner_caps_ = BuildCaps<Real>(cfg_.capOrder, cfg_.capNaz, cfg_.L, cfg_.R_in,  cfg_.y0, cfg_.z0, comm_);
    outer_caps_ = BuildCaps<Real>(cfg_.capOrder, cfg_.capNaz, cfg_.L, cfg_.R_out, cfg_.y0, cfg_.z0, comm_);
    inner_caps_.SetQuadScheme(quad_junctions::QJDefaultScheme<Real>(), cfg_.cov_q, cfg_.Nbeta, cfg_.max_depth);
    outer_caps_.SetQuadScheme(quad_junctions::QJDefaultScheme<Real>(), cfg_.cov_q, cfg_.Nbeta, cfg_.max_depth);

    // Concatenate node coords/normals in the SAME name-sorted order the operator
    // uses for its density blocks (see Assemble()): inner_caps, inner_tube,
    // outer_caps, outer_tube.
    sctl::Vector<Real> Xic, Xnic, Xit, Xnit, Xoc, Xnoc, Xot, Xnot;
    inner_caps_.GetNodeCoord(&Xic, &Xnic, nullptr);
    inner_tube_.GetNodeCoord(&Xit, &Xnit, nullptr);
    outer_caps_.GetNodeCoord(&Xoc, &Xnoc, nullptr);
    outer_tube_.GetNodeCoord(&Xot, &Xnot, nullptr);
    N_ic_ = Xic.Dim()/3; N_it_ = Xit.Dim()/3; N_oc_ = Xoc.Dim()/3; N_ot_ = Xot.Dim()/3;

    Xsurf_.ReInit(0); Xnsurf_.ReInit(0);
    auto app = [](sctl::Vector<Real>& dst, const sctl::Vector<Real>& src){ for (auto v : src) dst.PushBack(v); };
    app(Xsurf_, Xic); app(Xsurf_, Xit); app(Xsurf_, Xoc); app(Xsurf_, Xot);
    app(Xnsurf_, Xnic); app(Xnsurf_, Xnit); app(Xnsurf_, Xnoc); app(Xnsurf_, Xnot);

    // Per-node DL jump sign: -1 on the inner capsule (into-fluid normals),
    // +1 on the outer capsule (out-of-fluid normals).
    const sctl::Long Nn = N_ic_ + N_it_ + N_oc_ + N_ot_;
    jump_sign_.ReInit(Nn);
    for (sctl::Long i = 0;                 i < N_ic_ + N_it_; i++) jump_sign_[i] = (Real)-1;
    for (sctl::Long i = N_ic_ + N_it_;     i < Nn;            i++) jump_sign_[i] = (Real)+1;

    op_.reset();   // geometry changed -> drop any assembled operator (rebuilt lazily)
  }

  sctl::Long NumSurfNodes() const { return N_ic_ + N_it_ + N_oc_ + N_ot_; }
  const sctl::Vector<Real>& SurfCoord()  const { return Xsurf_; }
  const sctl::Vector<Real>& SurfNormal() const { return Xnsurf_; }
  const Config& GetConfig() const { return cfg_; }

  // Build the Dirichlet velocity BC (AoS, 3 per surface node): no-slip (0)
  // everywhere except the outer caps, where a Poiseuille-like axial (+x) profile
  // is prescribed -- inflow on the low-x dome, outflow on the high-x dome -- each
  // amplitude-normalized so its cross-cap volume flux magnitude equals |flux|.
  // On return, net_flux (should be ~0 = volume preservation) is reported.
  sctl::Vector<Real> BuildInflowOutflowBC(const Real flux, Real* net_flux_out = nullptr) const {
    // Amplitudes from the outer-cap far-field quadrature, split by axial half.
    // g = integral prof * (xhat . n) dA over each half; amp = flux / |g| gives
    // integral(v.n) = -flux (inflow, n.xhat<0) and +flux (outflow, n.xhat>0).
    Real g_in = 0, g_out = 0;
    {
      sctl::Vector<Real> X, Xn, wts, dist; sctl::Vector<sctl::Long> cnt;
      outer_caps_.GetFarFieldNodes(X, Xn, wts, dist, cnt, cfg_.tol);
      for (sctl::Long i = 0; i < wts.Dim(); i++) {
        const Real x = X[i*3+0];
        const Real prof = cap_profile(X[i*3+1], X[i*3+2]);
        const Real integrand = wts[i] * prof * Xn[i*3+0];   // prof * (xhat . n)
        if (x < cfg_.L/2) g_in += integrand; else g_out += integrand;
      }
    }
    g_in  = reduce_sum(g_in);
    g_out = reduce_sum(g_out);
    const Real amp_in  = (std::fabs(g_in)  > 0) ? flux / std::fabs(g_in)  : (Real)0;  // >0, +x inflow
    const Real amp_out = (std::fabs(g_out) > 0) ? flux / std::fabs(g_out) : (Real)0;  // >0, +x outflow

    const sctl::Long Nn = NumSurfNodes();
    sctl::Vector<Real> bc(Nn*3); bc.SetZero();
    // Outer-cap block occupies [N_ic_+N_it_+... ]? name-sorted order is
    // inner_caps, inner_tube, outer_caps, outer_tube -> outer_caps starts at
    // N_ic_+N_it_ and has N_oc_ nodes.
    const sctl::Long oc0 = N_ic_ + N_it_;
    for (sctl::Long n = oc0; n < oc0 + N_oc_; n++) {
      const Real x = Xsurf_[n*3+0], y = Xsurf_[n*3+1], z = Xsurf_[n*3+2];
      const Real prof = cap_profile(y, z);
      const Real amp = (x < cfg_.L/2) ? amp_in : amp_out;
      bc[n*3+0] = amp * prof;   // +x axial velocity; y,z components zero
    }

    // Net boundary flux integral(v.n) dA over the prescribed patches (far field).
    Real net = 0;
    {
      sctl::Vector<Real> X, Xn, wts, dist; sctl::Vector<sctl::Long> cnt;
      outer_caps_.GetFarFieldNodes(X, Xn, wts, dist, cnt, cfg_.tol);
      for (sctl::Long i = 0; i < wts.Dim(); i++) {
        const Real x = X[i*3+0];
        const Real amp = (x < cfg_.L/2) ? amp_in : amp_out;
        const Real prof = cap_profile(X[i*3+1], X[i*3+2]);
        net += wts[i] * amp * prof * Xn[i*3+0];   // v.n with v = amp*prof*xhat
      }
    }
    net = reduce_sum(net);
    if (net_flux_out) *net_flux_out = net;
    return bc;
  }

  // Build the moving-wall Dirichlet slip BC (AoS, 3 per surface node): on the
  // inner and outer TUBES, a purely radial slip velocity of magnitude rdot(x)
  // from the cached blended profiles (w(eta)*rdot_wall in the transition bands,
  // rdot_wall in the interior); no-slip (0) on the fixed dome caps. Mirrors
  // Annular::GetVslip: vslip = rdot/|X-Xc| * (X-Xc), with Xc = (x, y0, z0).
  // On return, net_flux_out reports  net = flux_outer - flux_inner  over the
  // tube surfaces; ~0 means the wall motion conserves the annular volume (the
  // interior-Dirichlet incompressibility compatibility condition). If it is not
  // ~0, add a compensating outer-cap flux (BuildInflowOutflowBC) equal to net.
  sctl::Vector<Real> BuildMovingWallSlipBC(Real* net_flux_out = nullptr) const {
    SCTL_ASSERT_MSG((bool)rdot_in_prof_ && (bool)rdot_out_prof_,
                    "BuildMovingWallSlipBC: no wall profile set (Config::r_*_of_x empty)");
    const sctl::Long Nn = NumSurfNodes();
    sctl::Vector<Real> bc(Nn*3); bc.SetZero();

    // Name-sorted node order is inner_caps, inner_tube, outer_caps, outer_tube.
    // Inner tube nodes: [N_ic_, N_ic_+N_it_). Outer tube: [N_ic_+N_it_+N_oc_, Nn).
    // The two cap blocks stay no-slip (0).
    const sctl::Long it0 = N_ic_;
    const sctl::Long ot0 = N_ic_ + N_it_ + N_oc_;
    auto set_slip = [&](sctl::Long n0, sctl::Long cnt, const std::function<Real(Real)>& rdot_of_x){
      for (sctl::Long n = n0; n < n0 + cnt; n++) {
        const Real x = Xsurf_[n*3+0], y = Xsurf_[n*3+1], z = Xsurf_[n*3+2];
        const Real ry = y - cfg_.y0, rz = z - cfg_.z0;
        const Real rn = std::sqrt(ry*ry + rz*rz);
        if (rn <= (Real)0) continue;
        const Real vr = rdot_of_x(x) / rn;   // rdot * unit-radial magnitude
        bc[n*3+0] = 0;                        // straight axis: no axial component
        bc[n*3+1] = vr * ry;
        bc[n*3+2] = vr * rz;
      }
    };
    set_slip(it0, N_it_, rdot_in_prof_);
    set_slip(ot0, N_ot_, rdot_out_prof_);

    if (net_flux_out)
      *net_flux_out = tube_flux(outer_tube_, rdot_out_prof_) - tube_flux(inner_tube_, rdot_in_prof_);
    return bc;
  }

  // Solve the interior combined-field BVP for the surface density sigma given the
  // Dirichlet velocity BC. SL_scal=-1, DL_scal=+1 (interior CFIE, no sign flip).
  sctl::Vector<Real> Solve(const sctl::Vector<Real>& bc, sctl::Long gmres_max_iter = 400,
                           const std::string& name = "capped-annulus") {
    StokesBIO<Real>& Op = Operator();   // assembled once; reused by EvaluateVelocity
    Op.SetTargetCoord(Xsurf_);          // on-surface targets for the self+jump matvec

    const Real jump = (Real)-0.5 * DL_scal_;   // interior
    const auto ApplyK = [&](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
      sctl::Vector<Real> Uc;
      Op.ComputePotential(Uc, sigma);          // SL_scal*S + DL_scal*D
      (*U) = Uc;
      for (sctl::Long n = 0; n < NumSurfNodes(); n++)
        for (int k = 0; k < 3; k++) (*U)[n*3+k] += jump_sign_[n] * jump * sigma[n*3+k];
    };

    sctl::GMRES<Real> solver(comm_, true);
    sctl::Vector<Real> sigma; sctl::Long iter = 0;
    const Real gmres_tol = cfg_.tol * (Real)10;
    solver(&sigma, ApplyK, bc, gmres_tol, gmres_max_iter, false, &iter);
    if (!comm_.Rank()) std::cout << "  " << name << ": GMRES iters=" << iter << "\n";

    // True relative residual.
    {
      sctl::Vector<Real> Ax; ApplyK(&Ax, sigma);
      Real rn = 0, bn = 0;
      for (sctl::Long i = 0; i < Ax.Dim(); i++) { const Real r = bc[i]-Ax[i]; rn += r*r; }
      for (sctl::Long i = 0; i < bc.Dim(); i++) bn += bc[i]*bc[i];
      rn = reduce_sum(rn);
      bn = reduce_sum(bn);
      if (!comm_.Rank())
        std::cout << "  " << name << ": TRUE rel residual = " << std::scientific
                  << std::sqrt(rn)/(bn>0?std::sqrt(bn):(Real)1) << std::defaultfloat << "\n";
    }
    return sigma;
  }

  // Evaluate the represented velocity field at off-surface interior targets Xtrg
  // (AoS). No jump term off-surface.
  void EvaluateVelocity(const sctl::Vector<Real>& sigma, const sctl::Vector<Real>& Xtrg,
                        sctl::Vector<Real>& U) {
    StokesBIO<Real>& Op = Operator();   // reuse the solve's assembled source quadrature
    Op.SetTargetCoord(Xtrg);            // off-surface targets: near-quad recomputed, self-quad reused
    Op.ComputePotential(U, sigma);
  }

  // Write the four surface geometries to VTU, tagging each file with its object
  // name. Each F_* is a per-node field on the corresponding block (AoS, 3 per
  // node) -- e.g. the split Dirichlet BC (see SplitSurfaceField). The caps are
  // static, so pass write_cap=true only once (typically step 0); moving tubes
  // are written every step. Note sctl's WriteVTK appends the ".vtu"/".pvtu"
  // extension itself, so we pass the prefix without one.
  void WriteVTK(const std::string& fname_prefix,
                const sctl::Vector<Real>& F_inner_tube, const sctl::Vector<Real>& F_outer_tube,
                const sctl::Vector<Real>& F_inner_caps, const sctl::Vector<Real>& F_outer_caps,
                const sctl::Comm& comm, bool write_cap = true) const {
    inner_tube_.WriteVTK(fname_prefix + "_inner_tube", F_inner_tube, comm);
    outer_tube_.WriteVTK(fname_prefix + "_outer_tube", F_outer_tube, comm);

    if (!write_cap) return;  // caps are static: skip unless requested
    inner_caps_.WriteVTK(fname_prefix + "_inner_caps", F_inner_caps, comm);
    outer_caps_.WriteVTK(fname_prefix + "_outer_caps", F_outer_caps, comm);
  }

  // Fill the four per-object VTUData blocks (inner tube, outer tube, inner caps,
  // outer caps), each colored by its per-node field F_* (AoS, 3 per node -- e.g.
  // the split Dirichlet BC from SplitSurfaceField). Unlike WriteVTK, the caps are
  // always emitted: a VTKHDF temporal step is self-contained, so each snapshot
  // must carry the full (static + moving) wall geometry. Used by the VTKHDF
  // time-series writer to merge all four surfaces into one grid per step.
  void GetSurfaceVTUData(sctl::VTUData& v_inner_tube, sctl::VTUData& v_outer_tube,
                         sctl::VTUData& v_inner_caps, sctl::VTUData& v_outer_caps,
                         const sctl::Vector<Real>& F_inner_tube, const sctl::Vector<Real>& F_outer_tube,
                         const sctl::Vector<Real>& F_inner_caps, const sctl::Vector<Real>& F_outer_caps) const {
    inner_tube_.GetVTUData(v_inner_tube, F_inner_tube);
    outer_tube_.GetVTUData(v_outer_tube, F_outer_tube);
    inner_caps_.GetVTUData(v_inner_caps, F_inner_caps);
    outer_caps_.GetVTUData(v_outer_caps, F_outer_caps);
  }

  // Split an AoS surface field defined over all NumSurfNodes() nodes (in the
  // fixed name-sorted order inner_caps, inner_tube, outer_caps, outer_tube -- the
  // same order the BC builders and the operator use) into the four per-geometry
  // blocks, in the argument order WriteVTK expects. Use this to turn a full
  // boundary-condition (or density/velocity) vector into per-object fields.
  void SplitSurfaceField(const sctl::Vector<Real>& F,
                         sctl::Vector<Real>& F_inner_tube, sctl::Vector<Real>& F_outer_tube,
                         sctl::Vector<Real>& F_inner_caps, sctl::Vector<Real>& F_outer_caps) const {
    SCTL_ASSERT_MSG(F.Dim() == NumSurfNodes()*3, "SplitSurfaceField: field size != 3*NumSurfNodes()");
    auto slice = [&](sctl::Long n0, sctl::Long cnt) {
      sctl::Vector<Real> out(cnt*3);
      for (sctl::Long i = 0; i < cnt*3; i++) out[i] = F[n0*3+i];
      return out;
    };
    F_inner_caps = slice(0,                     N_ic_);
    F_inner_tube = slice(N_ic_,                 N_it_);
    F_outer_caps = slice(N_ic_ + N_it_,         N_oc_);
    F_outer_tube = slice(N_ic_ + N_it_ + N_oc_, N_ot_);
  }

  // ParaView .pvd collection for the flow-geometry VTU series written by
  // WriteVTK, mirroring TracerParticles::WritePVD. Each timestep references that
  // step's moving inner/outer tube files plus the static cap files (written once
  // at step 0), so the caps stay visible throughout. `times` holds one timestamp
  // per written step; file names are step_prefix + zero-padded index + suffixes.
  // Static (geometry is recreated per step), so callers collect the times.
  static void WritePVD(const std::string& pvd_fname, const std::vector<Real>& times,
                       const std::string& step_prefix = "geom_step", const std::string& ext = ".pvtu") {
    std::ofstream pvd(pvd_fname);
    if (!pvd.is_open()) return;
    auto tag = [](sctl::Long i) { std::ostringstream os; os << std::setfill('0') << std::setw(4) << i; return os.str(); };

    pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";
    const std::string cap_in  = step_prefix + tag(0) + "_inner_caps" + ext;  // static caps (step 0)
    const std::string cap_out = step_prefix + tag(0) + "_outer_caps" + ext;
    for (sctl::Long i = 0; i < (sctl::Long)times.size(); ++i) {
      const std::string s = step_prefix + tag(i);
      pvd << "    <DataSet timestep=\"" << times[i] << "\" part=\"0\" file=\"" << s + "_inner_tube" + ext << "\"/>\n";
      pvd << "    <DataSet timestep=\"" << times[i] << "\" part=\"1\" file=\"" << s + "_outer_tube" + ext << "\"/>\n";
      pvd << "    <DataSet timestep=\"" << times[i] << "\" part=\"2\" file=\"" << cap_in  << "\"/>\n";
      pvd << "    <DataSet timestep=\"" << times[i] << "\" part=\"3\" file=\"" << cap_out << "\"/>\n";
    }
    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";
    pvd.close();
  }

 private:
  // Surface integral  ∮ rdot(x) dA  over a tube's far-field quadrature. On a
  // tube the moving-wall slip is purely radial, so u.n = rdot and this is the
  // tube's net volume-flux contribution.
  Real tube_flux(const sctl::SlenderElemList<Real>& tube, const std::function<Real(Real)>& rdot_of_x) const {
    sctl::Vector<Real> X, Xn, wts, dist; sctl::Vector<sctl::Long> cnt;
    tube.GetFarFieldNodes(X, Xn, wts, dist, cnt, cfg_.tol);
    Real s = 0;
    for (sctl::Long i = 0; i < wts.Dim(); i++) s += wts[i] * rdot_of_x(X[i*3+0]);
    return reduce_sum(s);
  }

  // Sum a local scalar across ranks (serial: identity).
  Real reduce_sum(Real x) const {
    sctl::StaticArray<Real,2> buf; buf[0] = x; buf[1] = 0;
    comm_.Allreduce((sctl::ConstIterator<Real>)buf, (sctl::Iterator<Real>)buf+1, 1, sctl::CommOp::SUM);
    return buf[1];
  }

  // Parabolic inflow/outflow profile on the OUTER cap, supported only on the
  // annular region R_in <= r <= R_out (r about the annulus centerline). The
  // r < R_in disk of the outer dome sits in the inner cap's shadow (blocked by
  // the inner capsule), so we prescribe zero there -- inflow/outflow is set only
  // on the outer tube cap's annular face, zero on the inner tube/cap. Same
  // profile drives both inflow (low-x) and outflow (high-x) caps.
  Real cap_profile(Real y, Real z) const {
    const Real dy = y - cfg_.y0, dz = z - cfg_.z0;
    const Real r2 = dy*dy + dz*dz;
    if (r2 < cfg_.R_in*cfg_.R_in) return (Real)0;   // inner-cap shadow: no flow
    const Real p = (Real)1 - r2 / (cfg_.R_out*cfg_.R_out);
    return p > 0 ? p : (Real)0;
  }

  // Add all four element lists in the fixed name-sorted order.
  void Assemble(StokesBIO<Real>& Op) const {
    Op.AddElemList(inner_caps_, "0_inner_caps");
    Op.AddElemList(inner_tube_, "1_inner_tube");
    Op.AddElemList(outer_caps_, "2_outer_caps");
    Op.AddElemList(outer_tube_, "3_outer_tube");
  }

  // Lazily assemble the combined-field operator once and reuse it for both the
  // GMRES matvec (on-surface targets) and off-surface velocity evaluation. The
  // (expensive) source self/near quadrature build thus happens a single time per
  // geometry; Setup() resets op_ so a moving-wall re-Setup rebuilds it.
  StokesBIO<Real>& Operator() {
    if (!op_) {
      op_ = std::make_shared<StokesBIO<Real>>(SL_scal_, DL_scal_, comm_);
      op_->SetAccuracy(cfg_.tol);
      Assemble(*op_);
    }
    return *op_;
  }

  Config cfg_;
  sctl::Comm comm_;
  const Real SL_scal_ = -1, DL_scal_ = +1;

  sctl::SlenderElemList<Real> inner_tube_, outer_tube_;
  sctl::QuadElemList<Real>    inner_caps_, outer_caps_;

  sctl::Vector<Real> Xsurf_, Xnsurf_, jump_sign_;
  sctl::Long N_ic_ = 0, N_it_ = 0, N_oc_ = 0, N_ot_ = 0;
  std::function<Real(Real)> rdot_in_prof_, rdot_out_prof_;  // blended rdot(x) for the slip BC
  std::shared_ptr<StokesBIO<Real>> op_;   // assembled-once combined-field operator
};

}  // namespace glymphbie

#endif  // _CAPPEDANNULUSFLOW_HPP_
