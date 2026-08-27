#ifndef _PERIODIC_BIO_HPP_
#define _PERIODIC_BIO_HPP_

// =============================================================================
// PeriodicBIO.hpp
//
// Support routines for the periodized combined-field Stokes BIE solve, matching
// the formulation in stokes-periodize-numtest (test/*.cpp). The periodic
// single-layer operator has a rank-one nullspace: the surface-mean of the
// density must be removed before applying the operator and added back to the
// resulting velocity. These helpers compute the (MPI-reduced) surface area and
// the surface-mean density over the inner + outer channel element lists.
//
//   u(x) = SL_scal * S[sigma0](x) + DL_scal * D[sigma0](x) + sigma_mean,
//   sigma0 = sigma - sigma_mean,
//   sigma_mean = (1/A) * \int_Gamma sigma dA.
// =============================================================================

#include "csbq/slender_element.hpp"
#include "csbq/slender_element.cpp"

// Quadrature-weighted surface integral of a dof-fast/node-slow field into I.
template <class Real> void SurfaceIntegral(sctl::Vector<Real>& I, const sctl::Vector<Real>& vals, const sctl::Vector<Real>& wts) {
  const sctl::Long dof = vals.Dim() / wts.Dim();
  SCTL_ASSERT(vals.Dim() == wts.Dim() * dof);
  if (I.Dim() != dof) I.ReInit(dof);
  I = 0;
  for (sctl::Long i = 0; i < wts.Dim(); i++) {
    for (sctl::Long j = 0; j < dof; j++) {
      I[j] += vals[i*dof + j] * wts[i];
    }
  }
}

// Add the per-component constant c0 to every node of a dof-fast/node-slow field.
template <class Real> void AddConstVec(sctl::Vector<Real>& vals, const sctl::Vector<Real>& c0) {
  const sctl::Long dof = c0.Dim();
  const sctl::Long N = vals.Dim() / dof;
  SCTL_ASSERT(vals.Dim() == N * dof);
  for (sctl::Long i = 0; i < N; i++) {
    for (sctl::Long j = 0; j < dof; j++) {
      vals[i*dof + j] += c0[j];
    }
  }
}

// Compute the far-field quadrature weights of the inner and outer channel
// element lists and the total surface area (summed across MPI processes).
template <class Real> void GetSurfWtsArea(const sctl::SlenderElemList<Real>& elem_inner,
                                          const sctl::SlenderElemList<Real>& elem_outer,
                                          sctl::Vector<Real>& wts_inner,
                                          sctl::Vector<Real>& wts_outer,
                                          Real& surface_area,
                                          const sctl::Comm& comm,
                                          const Real tol = 1e-14) {
  sctl::Vector<Real> X, Xn, dist_far;
  sctl::Vector<sctl::Long> element_wise_node_cnt;
  elem_inner.GetFarFieldNodes(X, Xn, wts_inner, dist_far, element_wise_node_cnt, tol);
  elem_outer.GetFarFieldNodes(X, Xn, wts_outer, dist_far, element_wise_node_cnt, tol);

  // Total (local) surface area = sum of far-field quadrature weights.
  sctl::Vector<Real> sa_loc(1); sa_loc[0] = 0;
  for (sctl::Long i = 0; i < wts_inner.Dim(); i++) sa_loc[0] += wts_inner[i];
  for (sctl::Long i = 0; i < wts_outer.Dim(); i++) sa_loc[0] += wts_outer[i];
  sctl::Vector<Real> sa_all(1); sa_all[0] = 0;
  comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
  surface_area = sa_all[0];
}

// Compute the surface-mean density sigma_mean (a 3-vector) over the inner +
// outer channel surfaces, reduced across MPI processes. <sigma> is the density
// on surface nodes ordered [inner; outer], with <size_inner> the number of
// inner-channel surface dofs. <wts_inner>/<wts_outer> and <surface_area> come
// from GetSurfWtsArea().
template <class Real> sctl::Vector<Real> ComputeSigmaMean(const sctl::Vector<Real>& sigma,
                                                          const sctl::SlenderElemList<Real>& elem_inner,
                                                          const sctl::SlenderElemList<Real>& elem_outer,
                                                          const sctl::Long size_inner,
                                                          const sctl::Vector<Real>& wts_inner,
                                                          const sctl::Vector<Real>& wts_outer,
                                                          const Real surface_area,
                                                          const sctl::Comm& comm) {
  // Split the density into inner/outer surface-node segments.
  sctl::Vector<Real> sigma_inner(size_inner, (sctl::Iterator<Real>) sigma.begin(), true);
  sctl::Vector<Real> sigma_outer(sigma.Dim() - size_inner, (sctl::Iterator<Real>) sigma.begin() + size_inner, true);

  // Interpolate to far-field quadrature nodes and integrate.
  sctl::Vector<Real> sigma_ff_in, sigma_ff_out;
  elem_inner.GetFarFieldDensity(sigma_ff_in, sigma_inner);
  elem_outer.GetFarFieldDensity(sigma_ff_out, sigma_outer);
  sctl::Vector<Real> I_in, I_out;
  SurfaceIntegral(I_in, sigma_ff_in, wts_inner);
  SurfaceIntegral(I_out, sigma_ff_out, wts_outer);

  sctl::Vector<Real> sa_loc(3);
  for (sctl::Long j = 0; j < 3; j++) sa_loc[j] = I_in[j] + I_out[j];

  // Reduce across MPI processes.
  sctl::Vector<Real> sigma_mean(3); sigma_mean = 0;
  for (sctl::Long j = 0; j < 3; j++) {
    comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin() + j, (sctl::Iterator<Real>) sigma_mean.begin() + j, 1, sctl::CommOp::SUM);
  }
  sigma_mean *= (1 / surface_area);
  return sigma_mean;
}

#endif
