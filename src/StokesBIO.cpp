#include "StokesBIO.hpp"

// StokesBIE by Dhairya Malhotra, 2025: 
// https://github.com/dmalhotra/stokes-periodize

template <class Real> 
StokesBIO<Real>::StokesBIO(const Real SL_scal, const Real DL_scal, const sctl::Comm comm)
  : comm_(comm), SL_scal_(SL_scal), DL_scal_(DL_scal), LayerPotenSL(ker_FxU, false, comm), LayerPotenDL(ker_DxU, false, comm) {
  LayerPotenSL.SetAccuracy(1e-14);
  LayerPotenDL.SetAccuracy(1e-14);
  LayerPotenSL.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  LayerPotenDL.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
};

template <class Real> 
void StokesBIO<Real>::SetPeriodicity(sctl::Periodicity periodicity, Real period_length) {
  LayerPotenSL.SetPeriodicity(periodicity, period_length);
  LayerPotenDL.SetPeriodicity(periodicity, period_length);
}

template <class Real> 
void StokesBIO<Real>::SetAccuracy(Real tol) {
  LayerPotenSL.SetAccuracy(tol);
  LayerPotenDL.SetAccuracy(tol);
}

template <class Real> 
void StokesBIO<Real>::AddElemList(const sctl::SlenderElemList<Real>& elem_lst, const std::string& name) {
  LayerPotenSL.AddElemList(elem_lst, name);
  LayerPotenDL.AddElemList(elem_lst, name);
}

template <class Real> 
void StokesBIO<Real>::DeleteElemList(const std::string& name) {
  LayerPotenSL.DeleteElemList(name);
  LayerPotenDL.DeleteElemList(name);
}

template <class Real> 
void StokesBIO<Real>::DeleteElemList() {
  LayerPotenSL.template DeleteElemList<sctl::SlenderElemList<Real>>();
  LayerPotenDL.template DeleteElemList<sctl::SlenderElemList<Real>>();
}

template <class Real> 
void StokesBIO<Real>::SetTargetCoord(const sctl::Vector<Real>& Xtrg) {
  LayerPotenSL.SetTargetCoord(Xtrg);
  LayerPotenDL.SetTargetCoord(Xtrg);
}

template <class Real> 
void StokesBIO<Real>::SetTargetNormal(const sctl::Vector<Real>& Xn_trg) {
  LayerPotenSL.SetTargetNormal(Xn_trg);
  LayerPotenDL.SetTargetNormal(Xn_trg);
}

template <class Real> 
sctl::Long StokesBIO<Real>::Dim(sctl::Integer k) const {
  return LayerPotenSL.Dim(k);
}

template <class Real> 
void StokesBIO<Real>::Setup() const {
  if (SL_scal_) LayerPotenSL.Setup();
  if (DL_scal_) LayerPotenDL.Setup();
}

template <class Real> 
void StokesBIO<Real>::ClearSetup() const {
  LayerPotenSL.ClearSetup();
  LayerPotenDL.ClearSetup();
}

template <class Real> 
void StokesBIO<Real>::ComputePotential(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) const {
  sctl::Vector<Real> Us, Ud;
  if (SL_scal_) LayerPotenSL.ComputePotential(Us, F);
  if (DL_scal_) LayerPotenDL.ComputePotential(Ud, F);

  if (SL_scal_ && DL_scal_) U = Us * SL_scal_ + Ud * DL_scal_;
  else if (SL_scal_) U = Us * SL_scal_;
  else if (DL_scal_) U = Ud * DL_scal_;
  else U.SetZero();
}

template <class Real> 
void StokesBIO<Real>::SqrtScaling(sctl::Vector<Real>& U) const {
  LayerPotenSL.SqrtScaling(U);
}

template <class Real> 
void StokesBIO<Real>::InvSqrtScaling(sctl::Vector<Real>& U) const {
  LayerPotenSL.InvSqrtScaling(U);
}

template class StokesBIO<float>;
template class StokesBIO<double>;



