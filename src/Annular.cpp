#include "Annular.hpp"

template <class Real>
Annular<Real>::Annular(
    const sctl::Vector<Real> Xc_inner_,
    const sctl::Vector<Real> Xc_outer_,
    const sctl::Vector<Real> r_inner_,
    const sctl::Vector<Real> r_outer_)
    : _Xc_inner(Xc_inner_),
      _Xc_outer(Xc_outer_),
      _r_inner(r_inner_),
      _r_outer(r_outer_)
{
    SCTL_ASSERT(_Xc_inner.Dim() == _r_inner.Dim() * 3);
    SCTL_ASSERT(_Xc_outer.Dim() == _r_outer.Dim() * 3);
}

template <class Real>
Annular<Real>::Annular(
    const sctl::Vector<Real> Xc_inner_,
    const sctl::Vector<Real> Xc_outer_,
    const sctl::Vector<Real> r_inner_,
    const sctl::Vector<Real> r_outer_,
    sctl::Comm comm_)
    : _Xc_inner(Xc_inner_),
      _Xc_outer(Xc_outer_),
      _r_inner(r_inner_),
      _r_outer(r_outer_),
      _comm(comm_)
{
    SCTL_ASSERT(_Xc_inner.Dim() == _r_inner.Dim() * 3);
    SCTL_ASSERT(_Xc_outer.Dim() == _r_outer.Dim() * 3);
}

template <class Real>
Annular<Real>::~Annular() {} 


template <class Real>
void Annular<Real>::SetInnerXc(const sctl::Vector<Real> X_in) {
    _Xc_inner = X_in;
    _SetupInner_bool = false;
}

template <class Real>
void Annular<Real>::SetOuterXc(const sctl::Vector<Real> X_in) {
    _Xc_outer = X_in;
    _SetupOuter_bool = false;
}

template <class Real>
void Annular<Real>::SetInnerR(const sctl::Vector<Real> r_in) {
    _r_inner = r_in;
    _SetupInner_bool = false;
}

template <class Real>
void Annular<Real>::SetOuterR(const sctl::Vector<Real> r_in) {
    _r_outer = r_in;
    _SetupOuter_bool = false;
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::Setup(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, 
    sctl::Vector<Real> drdt_outer_, bool force_setup) 
{
    if (force_setup) {
        _SetupInner_bool = false;
        _SetupOuter_bool = false;
    }
    if (_SetupInner_bool && _SetupOuter_bool) {
        return std::make_tuple(_Xc_inner, _Xc_outer, _r_inner, _r_outer, drdt_inner_, drdt_outer_);
    }
    sctl::Vector<Real> Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated, drdt_inner_updated, drdt_outer_updated;
    if (_SetupInner_bool) {
        std::tie(Xc_outer_updated, r_outer_updated, drdt_outer_updated) = SetupOuter(Nelem_, ElemOrder_, FourierOrder_, drdt_outer_);
        Xc_inner_updated = _Xc_inner;
        r_inner_updated = _r_inner;
        drdt_inner_updated = _drdt_inner;
    } else if (_SetupOuter_bool) {
        std::tie(Xc_inner_updated, r_inner_updated, drdt_inner_updated) = SetupInner(Nelem_, ElemOrder_, FourierOrder_, drdt_inner_);
        Xc_outer_updated = _Xc_outer;
        r_outer_updated = _r_outer;
        drdt_outer_updated = _drdt_outer;
    } else {
        std::tie(Xc_inner_updated, r_inner_updated, drdt_inner_updated) = SetupInner(Nelem_, ElemOrder_, FourierOrder_, drdt_inner_, force_setup);
        std::tie(Xc_outer_updated, r_outer_updated, drdt_outer_updated) = SetupOuter(Nelem_, ElemOrder_, FourierOrder_, drdt_outer_, force_setup);
    }    
    return std::make_tuple(Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated, drdt_inner_updated, drdt_outer_updated);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::Setup_mpi(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, 
    sctl::Vector<Real> drdt_outer_, bool force_setup) 
{
    if (_comm.Size() > Nelem_) {
        if (!_comm.Rank()) {
            std::cerr << "\n ERROR in Setup_mpi: Number of MPI processes larger than number of elements, will result in seg-faults. Use fewer processes." << std::endl;
        }
    }
    if (force_setup) {
        _SetupInner_bool = false;
        _SetupOuter_bool = false;
    }
    if (_SetupInner_bool && _SetupOuter_bool) {
        return std::make_tuple(_Xc_inner, _Xc_outer, _r_inner, _r_outer, drdt_inner_, drdt_outer_);
    }
    sctl::Vector<Real> Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated, drdt_inner_updated, drdt_outer_updated;
    if (_SetupInner_bool) {
        std::tie(Xc_outer_updated, r_outer_updated, drdt_outer_updated) = SetupOuter_mpi(Nelem_, ElemOrder_, FourierOrder_, drdt_outer_);
        Xc_inner_updated = _Xc_inner;
        r_inner_updated = _r_inner;
        drdt_inner_updated = _drdt_inner;
    } else if (_SetupOuter_bool) {
        std::tie(Xc_inner_updated, r_inner_updated, drdt_inner_updated) = SetupInner_mpi(Nelem_, ElemOrder_, FourierOrder_, drdt_inner_);
        Xc_outer_updated = _Xc_outer;
        r_outer_updated = _r_outer;
        drdt_outer_updated = _drdt_outer;
    } else {
        std::tie(Xc_inner_updated, r_inner_updated, drdt_inner_updated) = SetupInner_mpi(Nelem_, ElemOrder_, FourierOrder_, drdt_inner_, force_setup);
        std::tie(Xc_outer_updated, r_outer_updated, drdt_outer_updated) = SetupOuter_mpi(Nelem_, ElemOrder_, FourierOrder_, drdt_outer_, force_setup);
    }    
    return std::make_tuple(Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated, drdt_inner_updated, drdt_outer_updated);
}

template <class Real>
void Annular<Real>::InterpR(
    sctl::Vector<Real>& trg_r, const sctl::Vector<Real> src_r, 
    const sctl::Vector<Real> src_x, const sctl::Vector<Real> trg_x) 
{
    sctl::Matrix<Real> Minterp(src_x.Dim(), trg_x.Dim());
    sctl::Vector<Real> wts(src_x.Dim()*trg_x.Dim(), (sctl::Iterator<Real>)Minterp.begin(), false);
    sctl::LagrangeInterp<Real>::Interpolate(wts, src_x, trg_x); // row-major order, Ns x Nt (stacked) so wts for all targets from first source first, from second source, etc.
    
    // left multiply by source value of r to get row of target r values.
    SCTL_ASSERT(trg_r.Dim() == trg_x.Dim());
    sctl::Matrix<Real> trg_r_mat(1, trg_r.Dim(), (sctl::Iterator<Real>)trg_r.begin(), false); 
    sctl::Matrix<Real>::GEMM(trg_r_mat, sctl::Matrix<Real>(1, src_r.Dim(), (sctl::Iterator<Real>)src_r.begin(), false), Minterp);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::SetupInner(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, bool force_setup) 
{
    if (!_SetupInner_bool) {
        _Nelem_inner = Nelem_;
        _FourierOrder_inner = FourierOrder_;
        _ElemOrder_inner = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated, drdt_updated;
        if (!force_setup && _Xc_inner.Dim() == Nelem_ * ElemOrder_ * 3) {
            Xc_updated = _Xc_inner;
            r_updated = _r_inner;
            drdt_updated = drdt_inner_;
        } else {
            // Get (approximately) x nodes from original Xc
            sctl::Long OldOrder = _Xc_inner.Dim() / 3 / Nelem_; 
            sctl::Vector<Real> src_x(OldOrder);
            // Get new nodes
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            sctl::Vector<Real> trg_x(ElemOrder_);
            r_updated.ReInit(ElemOrder_ * Nelem_);
            drdt_updated.ReInit(ElemOrder_ * Nelem_);
            for (sctl::Long panel_ind=0; panel_ind < Nelem_; panel_ind ++) {
                src_x.SetZero();
                trg_x.SetZero();
                sctl::Vector<Real> Xc_updated_here(ElemOrder_*3, (sctl::Iterator<Real>) Xc_updated.begin() + panel_ind*ElemOrder_*3, false);
                sctl::Vector<Real> Xc_inner_here(OldOrder*3, (sctl::Iterator<Real>) _Xc_inner.begin() + panel_ind*OldOrder*3, false);
                for (sctl::Long node_ind = 0; node_ind < OldOrder; node_ind++) {
                    src_x[node_ind] = Xc_inner_here[node_ind*3];
                }
                for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
                    trg_x[node_ind] = Xc_updated_here[node_ind*3];
                }
                sctl::Vector<Real> r_inner_here(OldOrder, (sctl::Iterator<Real>) _r_inner.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> drdt_inner_here(OldOrder, (sctl::Iterator<Real>) drdt_inner_.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> r_updated_here(ElemOrder_, (sctl::Iterator<Real>) r_updated.begin() + panel_ind*ElemOrder_, false);
                sctl::Vector<Real> drdt_updated_here(ElemOrder_, (sctl::Iterator<Real>) drdt_updated.begin() + panel_ind*ElemOrder_, false);
                InterpR(r_updated_here, r_inner_here, src_x, trg_x);
                InterpR(drdt_updated_here, drdt_inner_here, src_x, trg_x);
            }
        }

        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_inner_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated);
        _elem_lst_inner = elem_lst_inner_;

        _SetupInner_bool = true;
        _Xc_inner = Xc_updated;
        _r_inner = r_updated;
        _drdt_inner = drdt_updated;
        _min_radius = _r_inner[0];
        for (const auto e : _r_inner) _min_radius = std::min<Real>(_min_radius, sctl::fabs(e));
        _elem_lst_inner.GetNodeCoord(&_X_inner, &_Xn_inner, nullptr);
    } else {
        std::cout << "Note: Inner already setup, nothing done." << std::endl;
    }
    return std::make_tuple(_Xc_inner, _r_inner, _drdt_inner);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::SetupInner_mpi(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, bool force_setup) 
{
    if (!_SetupInner_bool) {
        // Split ElemOrder up by number of mpi.
        // Assume all elements have the same ElemOrder and FourierOrder -- adaptive quadrature (later) should change only panel locations.
        sctl::Long Nmpi = _comm.Size();
        _loc_elem_cnt_inner = Nelem_ / Nmpi;
        sctl::Long remainder = Nelem_ - _loc_elem_cnt_inner * Nmpi;
        if (_comm.Rank() < remainder) {
            // if rank smaller than remainder from division, add a panel
            _loc_elem_cnt_inner += 1;
            _loc_elem_dsp_inner = _comm.Rank() * _loc_elem_cnt_inner;
        } else {
            _loc_elem_dsp_inner = _comm.Rank() * _loc_elem_cnt_inner + remainder;
        }
        _Nelem_inner = _loc_elem_cnt_inner;
        _FourierOrder_inner = FourierOrder_;
        _ElemOrder_inner = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated, drdt_updated;
        if (!force_setup && _Xc_inner.Dim() == _Nelem_inner * ElemOrder_ * 3) {
            Xc_updated = _Xc_inner;
            r_updated = _r_inner;
            drdt_updated = drdt_inner_;
        } else {
            // Get (approximately) x nodes from original Xc
            sctl::Long OldOrder = _Xc_inner.Dim() / 3 / _Nelem_inner; 
            sctl::Vector<Real> src_x(OldOrder);
            // Get new nodes
            std::tie(std::ignore, std::ignore, Xc_updated) = GetCenterLine_mpi(Nelem_, ElemOrder_, _comm); // Use overall Nelem_ in this function gives local Xc.
            SCTL_ASSERT(Xc_updated.Dim() == 3*_Nelem_inner*ElemOrder_);
            sctl::Vector<Real> trg_x(ElemOrder_);
            r_updated.ReInit(ElemOrder_ * _Nelem_inner);
            drdt_updated.ReInit(ElemOrder_ * _Nelem_inner);
            for (sctl::Long panel_ind=0; panel_ind < _Nelem_inner; panel_ind ++) {
                src_x.SetZero();
                trg_x.SetZero();
                sctl::Vector<Real> Xc_updated_here(ElemOrder_*3, (sctl::Iterator<Real>) Xc_updated.begin() + panel_ind*ElemOrder_*3, false);
                sctl::Vector<Real> Xc_inner_here(OldOrder*3, (sctl::Iterator<Real>) _Xc_inner.begin() + panel_ind*OldOrder*3, false);
                for (sctl::Long node_ind = 0; node_ind < OldOrder; node_ind++) {
                    src_x[node_ind] = Xc_inner_here[node_ind*3];
                }
                for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
                    trg_x[node_ind] = Xc_updated_here[node_ind*3];
                }
                sctl::Vector<Real> r_inner_here(OldOrder, (sctl::Iterator<Real>) _r_inner.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> drdt_inner_here(OldOrder, (sctl::Iterator<Real>) drdt_inner_.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> r_updated_here(ElemOrder_, (sctl::Iterator<Real>) r_updated.begin() + panel_ind*ElemOrder_, false);
                sctl::Vector<Real> drdt_updated_here(ElemOrder_, (sctl::Iterator<Real>) drdt_updated.begin() + panel_ind*ElemOrder_, false);
                InterpR(r_updated_here, r_inner_here, src_x, trg_x);
                InterpR(drdt_updated_here, drdt_inner_here, src_x, trg_x);
            }
        }

        // Again, assumes all elements have the same order.
        sctl::Vector<sctl::Long> ElemOrderVec(_loc_elem_cnt_inner);
        sctl::Vector<sctl::Long> FourierOrderVec(_loc_elem_cnt_inner);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_inner_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated);
        _elem_lst_inner = elem_lst_inner_;

        _SetupInner_bool = true;
        _Xc_inner = Xc_updated;
        _r_inner = r_updated;
        _drdt_inner = drdt_updated;
        _min_radius = _r_inner[0];
        for (const auto e : _r_inner) _min_radius = std::min<Real>(_min_radius, sctl::fabs(e));
        _elem_lst_inner.GetNodeCoord(&_X_inner, &_Xn_inner, nullptr);
    } else {
        std::cout << "Note: Inner already setup, nothing done." << std::endl;
    }
    return std::make_tuple(_Xc_inner, _r_inner, _drdt_inner);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::SetupOuter(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer_, bool force_setup) 
{
    if (!_SetupOuter_bool) {
        _Nelem_outer = Nelem_;
        _FourierOrder_outer = FourierOrder_;
        _ElemOrder_outer = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated, drdt_updated;
        if (!force_setup && _Xc_outer.Dim() == Nelem_ * ElemOrder_ * 3) {
            Xc_updated = _Xc_outer;
            r_updated = _r_outer;
            drdt_updated = drdt_outer_;
            // TODO: will add adaptive quadrature later. Also check what kind of adaptive quadrature is already implemented in CSBQ, should be some Fourier zooming in..?
        } else {
            // Get (approximately) x nodes from original Xc
            sctl::Long OldOrder = _Xc_outer.Dim() / 3 / Nelem_; 
            sctl::Vector<Real> src_x(OldOrder);
            // Get new nodes
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            sctl::Vector<Real> trg_x(ElemOrder_);
            r_updated.ReInit(ElemOrder_ * Nelem_);
            drdt_updated.ReInit(ElemOrder_ * Nelem_);
            for (sctl::Long panel_ind=0; panel_ind < Nelem_; panel_ind ++) {
                src_x.SetZero();
                trg_x.SetZero();
                sctl::Vector<Real> Xc_updated_here(ElemOrder_*3, (sctl::Iterator<Real>) Xc_updated.begin() + panel_ind*ElemOrder_*3, false);
                sctl::Vector<Real> Xc_outer_here(OldOrder*3, (sctl::Iterator<Real>) _Xc_outer.begin() + panel_ind*OldOrder*3, false);
                for (sctl::Long node_ind = 0; node_ind < OldOrder; node_ind++) {
                    src_x[node_ind] = Xc_outer_here[node_ind*3];
                }
                for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
                    trg_x[node_ind] = Xc_updated_here[node_ind*3];
                }
                sctl::Vector<Real> r_outer_here(OldOrder, (sctl::Iterator<Real>) _r_outer.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> drdt_outer_here(OldOrder, (sctl::Iterator<Real>) drdt_outer_.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> r_updated_here(ElemOrder_, (sctl::Iterator<Real>) r_updated.begin() + panel_ind*ElemOrder_, false);
                sctl::Vector<Real> drdt_updated_here(ElemOrder_, (sctl::Iterator<Real>) drdt_updated.begin() + panel_ind*ElemOrder_, false);
                InterpR(r_updated_here, r_outer_here, src_x, trg_x);
                InterpR(drdt_updated_here, drdt_outer_here, src_x, trg_x);
            }
        }
        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_outer_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated); // SlenderElemList will compute outward normal, so must be use Xn_outer to get correct orientation.
        _elem_lst_outer = elem_lst_outer_;

        _SetupOuter_bool = true;
        _Xc_outer = Xc_updated;
        _r_outer = r_updated;
        _drdt_outer = drdt_updated;
        _elem_lst_outer.GetNodeCoord(&_X_outer, &_Xn_outer, nullptr);
        // Correct normal of outer channel to point into fluid domain.
        _Xn_outer *= -1.;
    } else {
        std::cout << "Note: Outer already setup, nothing done." << std::endl;
    }
    return std::make_tuple(_Xc_outer, _r_outer, _drdt_outer);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
Annular<Real>::SetupOuter_mpi(
    const sctl::Long Nelem_, const sctl::Long ElemOrder_, 
    const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer_, bool force_setup) 
{
    if (!_SetupOuter_bool) {
        // Split ElemOrder up by number of mpi.
        // Assume all elements have the same ElemOrder and FourierOrder -- adaptive quadrature (later) should change only panel locations.
        sctl::Long Nmpi = _comm.Size();
        _loc_elem_cnt_outer = Nelem_ / Nmpi;
        sctl::Long remainder = Nelem_ - _loc_elem_cnt_outer * Nmpi;
        if (_comm.Rank() < remainder) {
            // if rank smaller than remainder from division, add a panel
            _loc_elem_cnt_outer += 1;
            _loc_elem_dsp_outer = _comm.Rank() * _loc_elem_cnt_outer;
        } else {
            _loc_elem_dsp_outer = _comm.Rank() * _loc_elem_cnt_outer + remainder;
        }

        _Nelem_outer = _loc_elem_cnt_outer;
        _FourierOrder_outer = FourierOrder_;
        _ElemOrder_outer = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated, drdt_updated;
        if (!force_setup && _Xc_outer.Dim() == _Nelem_outer * ElemOrder_ * 3) {
            Xc_updated = _Xc_outer;
            r_updated = _r_outer;
            drdt_updated = drdt_outer_;
            // TODO: will add adaptive quadrature later. Also check what kind of adaptive quadrature is already implemented in CSBQ, should be some Fourier zooming in..?
        } else {
            // Get (approximately) x nodes from original Xc
            sctl::Long OldOrder = _Xc_outer.Dim() / 3 / _Nelem_outer; 
            sctl::Vector<Real> src_x(OldOrder);
            // Get new nodes
            std::tie(std::ignore, std::ignore, Xc_updated) = GetCenterLine_mpi(Nelem_, ElemOrder_, _comm);
            SCTL_ASSERT(Xc_updated.Dim() == 3*_Nelem_outer*ElemOrder_);
            sctl::Vector<Real> trg_x(ElemOrder_);
            r_updated.ReInit(ElemOrder_ * _Nelem_outer);
            drdt_updated.ReInit(ElemOrder_ * _Nelem_outer);
            for (sctl::Long panel_ind=0; panel_ind < _Nelem_outer; panel_ind ++) {
                src_x.SetZero();
                trg_x.SetZero();
                sctl::Vector<Real> Xc_updated_here(ElemOrder_*3, (sctl::Iterator<Real>) Xc_updated.begin() + panel_ind*ElemOrder_*3, false);
                sctl::Vector<Real> Xc_outer_here(OldOrder*3, (sctl::Iterator<Real>) _Xc_outer.begin() + panel_ind*OldOrder*3, false);
                for (sctl::Long node_ind = 0; node_ind < OldOrder; node_ind++) {
                    src_x[node_ind] = Xc_outer_here[node_ind*3];
                }
                for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
                    trg_x[node_ind] = Xc_updated_here[node_ind*3];
                }
                sctl::Vector<Real> r_outer_here(OldOrder, (sctl::Iterator<Real>) _r_outer.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> drdt_outer_here(OldOrder, (sctl::Iterator<Real>) drdt_outer_.begin() + panel_ind*OldOrder, false);
                sctl::Vector<Real> r_updated_here(ElemOrder_, (sctl::Iterator<Real>) r_updated.begin() + panel_ind*ElemOrder_, false);
                sctl::Vector<Real> drdt_updated_here(ElemOrder_, (sctl::Iterator<Real>) drdt_updated.begin() + panel_ind*ElemOrder_, false);
                InterpR(r_updated_here, r_outer_here, src_x, trg_x);
                InterpR(drdt_updated_here, drdt_outer_here, src_x, trg_x);
            }
        }

        sctl::Vector<sctl::Long> ElemOrderVec(_loc_elem_cnt_outer);
        sctl::Vector<sctl::Long> FourierOrderVec(_loc_elem_cnt_outer);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_outer_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated); // SlenderElemList will compute outward normal, so must be use Xn_outer to get correct orientation.
        _elem_lst_outer = elem_lst_outer_;

        _SetupOuter_bool = true;
        _Xc_outer = Xc_updated;
        _r_outer = r_updated;
        _drdt_outer = drdt_updated;
        _elem_lst_outer.GetNodeCoord(&_X_outer, &_Xn_outer, nullptr);
        // Correct normal of outer channel to point into fluid domain.
        _Xn_outer *= -1.;
    } else {
        std::cout << "Note: Outer already setup, nothing done." << std::endl;
    }
    return std::make_tuple(_Xc_outer, _r_outer, _drdt_outer);
}

template <class Real>
void Annular<Real>::GetInnerCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out) {
    SCTL_ASSERT(_SetupInner_bool); 
    (*X_out) = _X_inner;
    if (Xn_out) {
        (*Xn_out) = _Xn_inner;
    }
}

template <class Real>
void Annular<Real>::GetOuterCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out) {
    SCTL_ASSERT(_SetupOuter_bool); 
    (*X_out) = _X_outer;
    if (Xn_out) {
        (*Xn_out) = _Xn_outer;
    }
}

template <class Real>
void Annular<Real>::GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn) {
    SCTL_ASSERT(_SetupInner_bool && _SetupOuter_bool); 
    if (Xn) {
        X->ReInit(_X_inner.Dim() + _X_outer.Dim());
        Xn->ReInit(_Xn_inner.Dim() + _Xn_outer.Dim());
        for (sctl::Long ind=0; ind<_X_inner.Dim(); ind++) {
            (*X)[ind] = _X_inner[ind];
            (*Xn)[ind] = _Xn_inner[ind];
        }
        sctl::Long offset = _X_inner.Dim();
        for (sctl::Long ind=0; ind<_X_outer.Dim(); ind++) {
            (*X)[ind+offset] = _X_outer[ind];
            (*Xn)[ind+offset] = _Xn_outer[ind];
        }
    } else {
        // no need to return Xn
        X->ReInit(_X_inner.Dim() + _X_outer.Dim());
        for (sctl::Long ind=0; ind<_X_inner.Dim(); ind++) {
            (*X)[ind] = _X_inner[ind];
        }
        sctl::Long offset = _X_inner.Dim();
        for (sctl::Long ind=0; ind<_X_outer.Dim(); ind++) {
            (*X)[ind+offset] = _X_outer[ind];
        }
    }
    
}

// x in [0,1], reinitializes x only, set y and z to 0.5... 
template <class Real>
sctl::Vector<Real> Annular<Real>::GetCenterLine(
    sctl::Long Nelem_, sctl::Long ElemOrder_) 
{
    // Returns centerline nodes when x=[0,1] divided into <Nelem> panels of <ElemOrder> number of Cheb nodes each
    const sctl::Vector<Real> nodes_01 = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder_);
    SCTL_ASSERT(ElemOrder_ == nodes_01.Dim());
    sctl::Vector<Real> Xc(ElemOrder_ * Nelem_ * 3);
    for (sctl::Long elem_ind = 0; elem_ind < Nelem_; elem_ind++) {
        for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
            const Real x = (elem_ind + nodes_01[node_ind]) / Nelem_; 
            Xc[(elem_ind*ElemOrder_+node_ind)*3 + 0] = x;
            Xc[(elem_ind*ElemOrder_+node_ind)*3 + 1] = 0.5;
            Xc[(elem_ind*ElemOrder_+node_ind)*3 + 2] = 0.5;
        }
    }
    return Xc;
}

template <class Real>
std::tuple<sctl::Long, sctl::Long, sctl::Vector<Real>> 
Annular<Real>::GetCenterLine_mpi(
    sctl::Long Nelem_, sctl::Long ElemOrder_, sctl::Comm comm_) 
{
    sctl::Long Nmpi = comm_.Size();
    sctl::Long loc_elem_cnt = Nelem_ / Nmpi;
    sctl::Long remainder = Nelem_ - loc_elem_cnt * Nmpi;
    sctl::Long loc_elem_dsp;
    if (comm_.Rank() < remainder) {
        // if rank smaller than remainder from division, add a panel
        loc_elem_cnt += 1;
        loc_elem_dsp = comm_.Rank() * loc_elem_cnt;
    } else {
        loc_elem_dsp = comm_.Rank() * loc_elem_cnt + remainder;
    }

    const sctl::Vector<Real> nodes_01 = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder_);
    sctl::Vector<Real> Xc(ElemOrder_ * loc_elem_cnt * 3);
    Xc = 0.5;
    for (sctl::Long elem_ind = loc_elem_dsp; elem_ind < loc_elem_dsp + loc_elem_cnt; elem_ind++) {
        for (sctl::Long node_ind = 0; node_ind < ElemOrder_; node_ind++) {
            const Real x = (elem_ind + nodes_01[node_ind]) / Nelem_; 
            Xc[((elem_ind - loc_elem_dsp)*ElemOrder_ + node_ind)*3 + 0] = x;
        }
    }

    return std::make_tuple(loc_elem_cnt, loc_elem_dsp, Xc);
}

template <class Real>
Real Annular<Real>::GetMinRadius() {
    SCTL_ASSERT(_SetupInner_bool);
    return _min_radius;
}

template <class Real>
sctl::Vector<sctl::Long> Annular<Real>::InDomain(sctl::Vector<Real> X_) {
    sctl::Long Nnodes = X_.Dim()/3;
    sctl::Vector<sctl::Long> in_domain(Nnodes);
    for (sctl::Long nd = 0; nd < Nnodes; nd++) {
        const Real x = X_[nd*3+0];
        const Real y = X_[nd*3+1];
        const Real z = X_[nd*3+2];
        if (x < _Xc_inner[0] || x < _Xc_outer[0] || x > _Xc_inner[_Xc_inner.Dim()-3] || x > _Xc_outer[_Xc_outer.Dim()-3]) {
            std::cout << "x value out of bounds , don't consider." << std::endl;
            in_domain[nd] = 0;
        } else {
            const auto search_panel = [](Real x, sctl::Long Nelem_, sctl::Long ElemOrder_, sctl::Vector<Real> Xc_) {
                for (sctl::Long ind=0; ind<Nelem_; ind++) {
                    if (x < Xc_[(ind+1)*ElemOrder_*3 - 3]) {
                        // x smaller than right end of this interval, return element index.
                        return ind; 
                    }
                }
                std::cout << "Did not find panel, something went wrong." << std::endl;
                return static_cast<sctl::Long>(-1);
            };
            // find panel x belongs to
            sctl::Long x_inner_ind = search_panel(x, _Nelem_inner, _ElemOrder_inner, _Xc_inner);
            sctl::Long x_outer_ind = search_panel(x, _Nelem_outer, _ElemOrder_outer, _Xc_outer);
            // interpolate r to this x value
            sctl::Vector<Real> r_inner_interp(1);
            sctl::Vector<Real> src_r(_ElemOrder_inner, (sctl::Iterator<Real>)_r_inner.begin()+x_inner_ind*_ElemOrder_inner, false);
            sctl::Vector<Real> src_x(_ElemOrder_inner);
            for (sctl::Long node_ind = 0; node_ind < _ElemOrder_inner; node_ind++) {
                src_x[node_ind] = _Xc_inner[x_inner_ind*_ElemOrder_inner*3 + node_ind*3];
            }
            sctl::Vector<Real> trg_x(1, (sctl::Iterator<Real>)X_.begin() + 3*nd, false);
            InterpR(r_inner_interp, src_r, src_x, trg_x);

            sctl::Vector<Real> r_outer_interp(1);
            src_x.SetZero();
            for (sctl::Long node_ind = 0; node_ind < _ElemOrder_outer; node_ind++) {
                src_x[node_ind] = _Xc_outer[x_outer_ind*_ElemOrder_outer*3 + node_ind*3];
            }
            sctl::Vector<Real> src_r2(_ElemOrder_outer, (sctl::Iterator<Real>)_r_outer.begin()+x_outer_ind*_ElemOrder_outer, false);
            InterpR(r_outer_interp, src_r2, src_x, trg_x);
            // compare r^2 with (y^2+z^2)
            in_domain[nd] = ((r_outer_interp[0]*r_outer_interp[0] - (y-0.5)*(y-0.5) - (z-0.5)*(z-0.5)) > 1e-5) && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5) - (r_inner_interp[0]*r_inner_interp[0]) > 1e-5);
        }
    }
    return in_domain;
}

template <class Real>
sctl::Vector<sctl::Long> Annular<Real>::InDomain_mpi(sctl::Vector<Real> X_) {
    sctl::Long Nnodes = X_.Dim()/3;
    sctl::Vector<sctl::Long> in_domain(Nnodes);
    for (sctl::Long nd = 0; nd < Nnodes; nd++) {
        const Real x = X_[nd*3+0];
        const Real y = X_[nd*3+1];
        const Real z = X_[nd*3+2];
        if (x < _X_inner[0] || x < _X_outer[0] || x > _X_inner[_X_inner.Dim()-3] || x > _X_outer[_X_outer.Dim()-3]) {
            std::cout << "x value out of bounds for process rank " << _comm.Rank() << ", don't consider." << std::endl;
            in_domain[nd] = 0;
        } else {
            const auto search_panel = [](Real x, sctl::Long loc_elem_dsp_, sctl::Long loc_elem_cnt_, sctl::Long ElemOrder_, sctl::Vector<Real> Xc_) {
                for (sctl::Long ind=0; ind<loc_elem_cnt_; ind++) {
                    if (x < Xc_[(ind+1)*ElemOrder_*3 - 3]) {// extra assumption that all elements have order ElemOrder_
                        // x smaller than right end of this interval, return  element index.
                        return ind; 
                    }
                }
                std::cout << "Did not find panel, something went wrong." << std::endl;
                return static_cast<sctl::Long>(-1);
            };
            // find panel x belongs to
            sctl::Long x_inner_ind = search_panel(x, _loc_elem_dsp_inner, _loc_elem_cnt_inner, _ElemOrder_inner, _Xc_inner);
            sctl::Long x_outer_ind = search_panel(x, _loc_elem_dsp_outer, _loc_elem_cnt_outer, _ElemOrder_outer, _Xc_outer);
            // interpolate r to this x value
            sctl::Vector<Real> r_inner_interp(1);
            sctl::Vector<Real> src_r(_ElemOrder_inner, (sctl::Iterator<Real>)_r_inner.begin()+x_inner_ind*_ElemOrder_inner, false);
            sctl::Vector<Real> src_x(_ElemOrder_inner);
            for (sctl::Long node_ind = 0; node_ind < _ElemOrder_inner; node_ind++) {
                src_x[node_ind] = _Xc_inner[x_inner_ind*_ElemOrder_inner*3 + node_ind*3];
            }
            sctl::Vector<Real> trg_x(1, (sctl::Iterator<Real>)X_.begin() + 3*nd, false);
            InterpR(r_inner_interp, src_r, src_x, trg_x);

            sctl::Vector<Real> r_outer_interp(1);
            src_x.SetZero();
            for (sctl::Long node_ind = 0; node_ind < _ElemOrder_outer; node_ind++) {
                src_x[node_ind] = _Xc_outer[x_outer_ind*_ElemOrder_outer*3 + node_ind*3];
            }
            sctl::Vector<Real> src_r2(_ElemOrder_outer, (sctl::Iterator<Real>)_r_outer.begin()+x_outer_ind*_ElemOrder_outer, false);
            InterpR(r_outer_interp, src_r2, src_x, trg_x);
            // compare r^2 with (y^2+z^2)
            in_domain[nd] = ((r_outer_interp[0]*r_outer_interp[0] - (y-0.5)*(y-0.5) - (z-0.5)*(z-0.5)) > 1e-5) && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5) - (r_inner_interp[0]*r_inner_interp[0]) > 1e-5);
        }
    }
    return in_domain;
}

template <class Real>
sctl::SlenderElemList<Real> Annular<Real>::GetInnerElemList() {
    return _elem_lst_inner;
}

template <class Real>
sctl::SlenderElemList<Real> Annular<Real>::GetOuterElemList() {
    return _elem_lst_outer;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslip() {
    // Given drdt at the same x discre as centerline nodes, return the axisymmetric vslip at surface nodes.
    sctl::Vector<Real> v1, v2;
    v1 = GetVslipInner();
    v2 = GetVslipOuter();
    sctl::Vector<Real> v(v1.Dim() + v2.Dim());
    for (sctl::Long ind=0; ind<v.Dim(); ind++) {
        if (ind < v1.Dim()) {
            v[ind] = v1[ind];
        } else {
            sctl::Long ind_v2 = ind - v1.Dim();
            v[ind] = v2[ind_v2];
        }
    }
    return v;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslip_mpi() {
    // Given drdt at the same x discre as centerline nodes, return the axisymmetric vslip at surface nodes.
    sctl::Vector<Real> v1, v2;
    v1 = GetVslipInner_mpi();
    v2 = GetVslipOuter_mpi();
    sctl::Vector<Real> v(v1.Dim() + v2.Dim());
    for (sctl::Long ind=0; ind<v.Dim(); ind++) {
        if (ind < v1.Dim()) {
            v[ind] = v1[ind];
        } else {
            sctl::Long ind_v2 = ind - v1.Dim();
            v[ind] = v2[ind_v2];
        }
    }
    return v;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslipInner() {
    SCTL_ASSERT(_SetupInner_bool); // check that inner channel is set up.
    // r vector for each angular node, per center node.
    sctl::Vector<Real> vslip(_X_inner.Dim());
    SCTL_ASSERT(vslip.Dim() == _Nelem_inner*_ElemOrder_inner*_FourierOrder_inner * 3);

    for (sctl::Long cn=0; cn < _Nelem_inner * _ElemOrder_inner; cn++) {
        sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) _Xc_inner.begin() + cn*3, false); // iterator starting at the x coord of current center node.
        for (sctl::Long fn = 0; fn < _FourierOrder_inner; fn++) { // check that this is the number of F nodes
            sctl::Long fncn = 3*(cn * _FourierOrder_inner + fn);
            sctl::Vector<Real> X_here(3, (sctl::Iterator<Real>) _X_inner.begin() + fncn, false); 
            sctl::Vector<Real> rvec_here = X_here - Xc_here;
            Real rnorm = rvec_here[0]*rvec_here[0] + rvec_here[1]*rvec_here[1] + rvec_here[2]*rvec_here[2];
            rnorm = sctl::sqrt<Real>(rnorm);
            sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) vslip.begin() + fncn, false);
            vslip_here = _drdt_inner[cn] / rnorm * rvec_here; // vslip points out in r direction with magnitude drdt.
        }
    }
    return vslip;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslipInner_mpi() {
    SCTL_ASSERT(_SetupInner_bool); // check that inner channel is set up.
    // r vector for each angular node, per center node.
    sctl::Vector<Real> vslip(_X_inner.Dim());
    SCTL_ASSERT(vslip.Dim() == _loc_elem_cnt_inner*_ElemOrder_inner*_FourierOrder_inner * 3);

    for (sctl::Long cn=0; cn < _loc_elem_cnt_inner * _ElemOrder_inner; cn++) {
        sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) _Xc_inner.begin() + cn*3, false); // iterator starting at the x coord of current center node.
        for (sctl::Long fn = 0; fn < _FourierOrder_inner; fn++) { // check that this is the number of F nodes
            sctl::Long fncn = 3*(cn * _FourierOrder_inner + fn);
            sctl::Vector<Real> X_here(3, (sctl::Iterator<Real>) _X_inner.begin() + fncn, false); 
            sctl::Vector<Real> rvec_here = X_here - Xc_here;
            Real rnorm = rvec_here[0]*rvec_here[0] + rvec_here[1]*rvec_here[1] + rvec_here[2]*rvec_here[2];
            rnorm = sctl::sqrt<Real>(rnorm);
            sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) vslip.begin() + fncn, false);
            vslip_here = _drdt_inner[cn] / rnorm * rvec_here; // vslip points out in r direction with magnitude drdt.
        }
    }
    return vslip;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslipOuter() {
    SCTL_ASSERT(_SetupOuter_bool); // check that outer channel is set up.
    // r vector for each angular node, per center node.
    sctl::Vector<Real> vslip(_X_outer.Dim());
    SCTL_ASSERT(vslip.Dim() == _Nelem_outer*_ElemOrder_outer*_FourierOrder_outer * 3);

    for (sctl::Long cn=0; cn < _Nelem_outer * _ElemOrder_outer; cn++) {
        sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) _Xc_outer.begin() + cn*3, false); // iterator starting at the x coord of current center node.
        for (sctl::Long fn = 0; fn < _FourierOrder_outer; fn++) { // check that this is the number of F nodes
            sctl::Long fncn = 3*(cn * _FourierOrder_outer + fn);
            sctl::Vector<Real> X_here(3, (sctl::Iterator<Real>) _X_outer.begin() + fncn, false); 
            sctl::Vector<Real> rvec_here = X_here - Xc_here;
            Real rnorm = rvec_here[0]*rvec_here[0] + rvec_here[1]*rvec_here[1] + rvec_here[2]*rvec_here[2];
            rnorm = sctl::sqrt<Real>(rnorm);
            sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) vslip.begin() + fncn, false);
            vslip_here = _drdt_outer[cn] / rnorm * rvec_here; // vslip points out in r direction with magnitude drdt.
        }
    }
    return vslip;
}

template <class Real>
sctl::Vector<Real> Annular<Real>::GetVslipOuter_mpi() {
    SCTL_ASSERT(_SetupOuter_bool); // check that outer channel is set up.
    // r vector for each angular node, per center node.
    sctl::Vector<Real> vslip(_X_outer.Dim());
    SCTL_ASSERT(vslip.Dim() == _loc_elem_cnt_outer*_ElemOrder_outer*_FourierOrder_outer * 3);

    for (sctl::Long cn=0; cn < _loc_elem_cnt_outer * _ElemOrder_outer; cn++) {
        sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) _Xc_outer.begin() + cn*3, false); // iterator starting at the x coord of current center node.
        for (sctl::Long fn = 0; fn < _FourierOrder_outer; fn++) { // check that this is the number of F nodes
            sctl::Long fncn = 3*(cn * _FourierOrder_outer + fn);
            sctl::Vector<Real> X_here(3, (sctl::Iterator<Real>) _X_outer.begin() + fncn, false); 
            sctl::Vector<Real> rvec_here = X_here - Xc_here;
            Real rnorm = rvec_here[0]*rvec_here[0] + rvec_here[1]*rvec_here[1] + rvec_here[2]*rvec_here[2];
            rnorm = sctl::sqrt<Real>(rnorm);
            sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) vslip.begin() + fncn, false);
            vslip_here = _drdt_outer[cn] / rnorm * rvec_here; // vslip points out in r direction with magnitude drdt.
        }
    }
    return vslip;
}

template <class Real>
void Annular<Real>::WriteVTK(std::string filename_inner, std::string filename_outer, sctl::Vector<Real> F_inner, sctl::Vector<Real> F_outer, sctl::Comm comm_) {
    _elem_lst_inner.WriteVTK(filename_inner, F_inner, comm_);
    _elem_lst_outer.WriteVTK(filename_outer, F_outer, comm_);
}

template class Annular<float>;
template class Annular<double>;
