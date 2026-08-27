#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

/**
 * Particle.hpp
 *
 * Lagrangian tracer particles carried by the annular perivascular flow.
 *
 * TracerParticles<Real> holds a set of massless tracers by position (AoS
 * {x1,y1,z1,x2,...}), records their full position history as a time series, and
 * advects them through a caller-supplied velocity field. It knows the annular
 * configuration (axial extent [x0,x1]=[0,L], centerline (y0,z0), inner/outer
 * radii R_in,R_out) so it can:
 *   - test whether a tracer is inside the annular fluid gap (InDomain),
 *   - advance positions with an explicit step given a velocity evaluator, and
 *   - apply the periodic re-seed at the outlet: a tracer that passes the outlet
 *     cap (x > x1) is copied back to the inlet at the SAME radial/azimuthal
 *     position (x -> x - L, keeping y,z), treating the finite capped annulus as
 *     one period of an infinite perivascular space. If the moved walls make the
 *     re-seeded point non-fluid, its radius is clamped back into (R_in,R_out).
 *
 * Output: WriteVTK writes one VTU point cloud per timestep (VTK_VERTEX cells,
 * colored by velocity and particle id) for flow snapshots; WriteTrajectoryVTK
 * writes the whole recorded history as one VTK_POLY_LINE per particle.
 *
 * Header-only; template on Real. The velocity evaluator is any callable
 *   void(const sctl::Vector<Real>& X, sctl::Vector<Real>& U)
 * (e.g. a lambda wrapping CappedAnnulusFlow::EvaluateVelocity), so this class is
 * independent of the specific BIE solver.
 */

#include <sctl.hpp>
#include <vector>
#include <string>
#include <cmath>

template <class Real>
class TracerParticles {
 public:
  struct Geometry {
    Real x0 = 0.0, x1 = 1.0;     // axial extent [inlet, outlet]
    Real y0 = 0.5, z0 = 0.5;     // annulus centerline (y,z)
    Real R_in = 0.3, R_out = 0.48;
    Real in_tol = 1e-5;          // radial margin for the in-domain test
  };

  TracerParticles() = default;
  TracerParticles(const sctl::Vector<Real>& X0, const Geometry& geom)
      : _X(X0), _geom(geom) {}

  // Replace the tracer positions (AoS). Clears nothing else.
  void SetPositions(const sctl::Vector<Real>& X) { _X = X; }
  const sctl::Vector<Real>& Positions() const { return _X; }
  sctl::Long Size() const { return _X.Dim()/3; }
  void SetGeometry(const Geometry& g) { _geom = g; }
  const Geometry& GetGeometry() const { return _geom; }

  // Per-particle in-domain mask: 1 if inside the annular fluid gap, else 0.
  sctl::Vector<sctl::Long> InDomain() const { return InDomain(_X); }
  sctl::Vector<sctl::Long> InDomain(const sctl::Vector<Real>& X) const {
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<sctl::Long> in(N);
    for (sctl::Long i = 0; i < N; i++) {
      const Real x = X[i*3+0], dy = X[i*3+1]-_geom.y0, dz = X[i*3+2]-_geom.z0;
      const Real r2 = dy*dy + dz*dz;
      const bool radial = (r2 < _geom.R_out*_geom.R_out - _geom.in_tol) &&
                          (r2 > _geom.R_in*_geom.R_in  + _geom.in_tol);
      const bool axial  = (x >= _geom.x0) && (x <= _geom.x1);
      in[i] = (radial && axial) ? 1 : 0;
    }
    return in;
  }

  // Advance one explicit step X += dt*U using the supplied velocity evaluator,
  // then apply the outlet->inlet periodic re-seed. VelFn signature:
  //   void(const sctl::Vector<Real>& X, sctl::Vector<Real>& U)
  template <class VelFn>
  void Step(Real dt, VelFn&& vel) {
    sctl::Vector<Real> U;
    vel(_X, U);
    SCTL_ASSERT_MSG(U.Dim() == _X.Dim(), "TracerParticles::Step velocity/position size mismatch");
    for (sctl::Long i = 0; i < _X.Dim(); i++) _X[i] += dt * U[i];
    ReseedAtInlet();
  }

  // Outlet -> inlet periodic transition: a tracer past the outlet cap (x > x1)
  // is copied to the inlet (x -= L) at the same (y,z); if the result is not in
  // the fluid gap, its radius is clamped into (R_in, R_out).
  void ReseedAtInlet() {
    const Real L = _geom.x1 - _geom.x0;
    for (sctl::Long i = 0; i < Size(); i++) {
      Real& x = _X[i*3+0];
      if (x > _geom.x1) {
        x -= L;
        ClampRadius(i);
      } else if (x < _geom.x0) {   // symmetric: re-enter from the outlet end
        x += L;
        ClampRadius(i);
      }
    }
  }

  // Append the current positions to the history with timestamp t.
  void RecordHistory(Real t) {
    _X_hist.push_back(_X);
    _times.push_back(t);
  }
  sctl::Long NumSnapshots() const { return (sctl::Long)_X_hist.size(); }

  // Build the current positions as a VTU point cloud (VTK_VERTEX), colored by
  // the velocity U (AoS, optional). Shared by WriteVTK and the VTKHDF writer.
  void GetVTUData(sctl::VTUData& vtu, const sctl::Vector<Real>& U = sctl::Vector<Real>()) const {
    const sctl::Long N = Size();
    for (sctl::Long i = 0; i < N; i++)
      for (sctl::Integer k = 0; k < 3; k++) vtu.coord.PushBack((sctl::VTUData::VTKReal)_X[i*3+k]);
    // velocity (zeros if not provided) as a 3-vector value
    const bool haveU = (U.Dim() == N*3);
    for (sctl::Long i = 0; i < N; i++)
      for (sctl::Integer k = 0; k < 3; k++)
        vtu.value.PushBack((sctl::VTUData::VTKReal)(haveU ? U[i*3+k] : (Real)0));
    for (sctl::Long i = 0; i < N; i++) {
      vtu.connect.PushBack((int32_t)i);
      vtu.offset.PushBack((int32_t)(i+1));
      vtu.types.PushBack((uint8_t)1);   // VTK_VERTEX
    }
  }

  // Write the current positions as a VTU point cloud (VTK_VERTEX), colored by
  // the velocity U (AoS, optional) and particle id. `fname` gets ".vtu".
  void WriteVTK(const std::string& fname, const sctl::Vector<Real>& U = sctl::Vector<Real>()) const {
    sctl::VTUData vtu;
    GetVTUData(vtu, U);
    vtu.WriteVTK(fname, sctl::Comm::Self());
  }

  // Write the recorded history as one polyline per particle (VTK_POLY_LINE),
  // giving each tracer's trajectory in time. Requires >= 1 recorded snapshot.
  void WriteTrajectoryVTK(const std::string& fname) const {
    sctl::VTUData vtu;
    const sctl::Long Ns = NumSnapshots();
    if (Ns == 0) return;
    const sctl::Long N = _X_hist[0].Dim()/3;   // particle count at first snapshot
    // Points ordered [particle][snapshot]; time stored as the scalar value.
    for (sctl::Long i = 0; i < N; i++)
      for (sctl::Long s = 0; s < Ns; s++)
        for (sctl::Integer k = 0; k < 3; k++)
          vtu.coord.PushBack((sctl::VTUData::VTKReal)_X_hist[s][i*3+k]);
    for (sctl::Long i = 0; i < N; i++)
      for (sctl::Long s = 0; s < Ns; s++)
        vtu.value.PushBack((sctl::VTUData::VTKReal)_times[s]);
    for (sctl::Long i = 0; i < N; i++) {
      for (sctl::Long s = 0; s < Ns; s++) vtu.connect.PushBack((int32_t)(i*Ns + s));
      vtu.offset.PushBack((int32_t)((i+1)*Ns));
      vtu.types.PushBack((uint8_t)4);   // VTK_POLY_LINE
    }
    vtu.WriteVTK(fname, sctl::Comm::Self());
  }

  // Generates a ParaView .pvd collection file referencing all recorded snapshots
  void WritePVD(const std::string& pvd_fname, 
                const std::string& step_prefix = "tracers_step", 
                const std::string& ext = ".vtu") const {
    std::ofstream pvd(pvd_fname);
    if (!pvd.is_open()) return;

    pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";

    for (sctl::Long i = 0; i < NumSnapshots(); ++i) {
      std::ostringstream ss;
      // Formats step index into 4 zero-padded digits (e.g., tracers_step0001.vtu)
      ss << step_prefix << std::setfill('0') << std::setw(4) << i << ext;
      
      pvd << "    <DataSet timestep=\"" << _times[i] 
          << "\" file=\"" << ss.str() << "\"/>\n";
    }

    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";
    pvd.close();
  }

 private:
  // Clamp particle i's radius into the open gap (R_in, R_out), keeping its angle.
  void ClampRadius(sctl::Long i) {
    const Real dy = _X[i*3+1]-_geom.y0, dz = _X[i*3+2]-_geom.z0;
    const Real r = std::sqrt(dy*dy + dz*dz);
    const Real lo = _geom.R_in + 10*_geom.in_tol, hi = _geom.R_out - 10*_geom.in_tol;
    if (r <= 0) { _X[i*3+1] = _geom.y0 + (Real)0.5*(lo+hi); _X[i*3+2] = _geom.z0; return; }
    Real rc = r; if (rc < lo) rc = lo; if (rc > hi) rc = hi;
    if (rc != r) { const Real s = rc/r; _X[i*3+1] = _geom.y0 + dy*s; _X[i*3+2] = _geom.z0 + dz*s; }
  }

  sctl::Vector<Real> _X;                     // current positions, AoS
  std::vector<sctl::Vector<Real>> _X_hist;   // position history (time series)
  std::vector<Real> _times;                  // timestamps of each snapshot
  Geometry _geom;
};

#endif  // _PARTICLE_HPP_
