#include "Annular.hpp"

// TODO: Implement OPENMP AND MPI (AND MAYBE KOKKOS) parallelization.
// TODO: In setup, interpolating for radius at new centerline nodes assume r = r(x), no y and z dependency.
// NOTE: Reinitialized centerline is at (x, 0.5, 0.5).

template <class Real>
Annular<Real>::Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_) : \
Xc_inner(Xc_inner_), Xc_outer(Xc_outer_), r_inner(r_inner_), r_outer(r_outer_)
{
    SCTL_ASSERT(Xc_inner_.Dim()==r_inner_.Dim()*3);
    SCTL_ASSERT(Xc_outer_.Dim()==r_outer_.Dim()*3);
}

template <class Real>
Annular<Real>::~Annular() {} 


template <class Real>
void Annular<Real>::SetInnerXc(const sctl::Vector<Real> X_in) {
    Xc_inner = X_in;
    SetupInner_bool = false;
}

template <class Real>
void Annular<Real>::SetOuterXc(const sctl::Vector<Real> X_in) {
    Xc_outer = X_in;
    SetupOuter_bool = false;
}

template <class Real>
void Annular<Real>::SetInnerR(const sctl::Vector<Real> r_in) {
    r_inner = r_in;
    SetupInner_bool = false;
}

template <class Real>
void Annular<Real>::SetOuterR(const sctl::Vector<Real> r_in) {
    r_outer = r_in;
    SetupOuter_bool = false;
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Annular<Real>::Setup(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_) {
    if (SetupInner_bool && SetupOuter_bool) {
        // Have already Setup() before, only redo setup if Nelem, ElemOrder, or FourierOrder are different
        if (Nelem_==Nelem_inner && Nelem_==Nelem_outer && ElemOrder_==ElemOrder_inner && ElemOrder_==ElemOrder_outer && FourierOrder_==FourierOrder_inner && FourierOrder_== FourierOrder_outer) {
            std::cout << "Already set up with the same parameters, nothing done." << std::endl;
            return std::make_tuple(Xc_inner, Xc_outer, r_inner, r_outer);
        }
    }
    sctl::Vector<Real> Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated;
    if (SetupInner_bool) {
        // Inner set up, outer not
        std::cout << "Inner is setup, setting up outer now."<<std::endl;
        auto [Xc_outer_updated, r_outer_updated] = SetupOuter(Nelem_, ElemOrder_, FourierOrder_);
        Xc_inner_updated = Xc_inner;
        r_inner_updated = r_inner;
    } else if (SetupOuter_bool) {
        // Outer set up, inner not
        std::cout << "OUter is setup, setting up inner now." << std::endl;
        auto [Xc_inner_updated, r_inner_updated] = SetupInner(Nelem_, ElemOrder_, FourierOrder_);
        Xc_outer_updated = Xc_outer;
        r_outer_updated = r_outer;
    } else {
        // Neither is set up
        std::cout << "neither is setup." << std::endl;
        auto [Xc_inner_updated, r_inner_updated] = SetupInner(Nelem_, ElemOrder_, FourierOrder_);
        auto [Xc_outer_updated, r_outer_updated] = SetupOuter(Nelem_, ElemOrder_, FourierOrder_);
    }    
    return std::make_tuple(Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated);
}

template <class Real>
void Annular<Real>::InterpR(sctl::Vector<Real>& trg_r, const sctl::Vector<Real> src_r, const sctl::Vector<Real> src_x, const sctl::Vector<Real> trg_x) {
    sctl::Matrix<Real> Minterp(src_x.Dim(), trg_x.Dim());
    sctl::Vector<Real> wts(src_x.Dim()*trg_x.Dim(), (sctl::Iterator<Real>)Minterp.begin(), false);
    sctl::LagrangeInterp<Real>::Interpolate(wts, src_x, trg_x); // row-major order, Ns x Nt (stacked) so wts for all targets from first source first, from second source, etc.
    // left multiply by source value of r to get row of target r values.
    trg_r.ReInit(trg_x.Dim());
    sctl::Matrix<Real> trg_r_mat(1, trg_r.Dim(), (sctl::Iterator<Real>)trg_r.begin(),false); 
    sctl::Matrix<Real>::GEMM(trg_r_mat, sctl::Matrix<Real>(1,r_inner.Dim(),(sctl::Iterator<Real>)r_inner.begin(),false), Minterp);
    // TODO: check that this assignes values to trg_r
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> Annular<Real>::SetupInner(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_) {
    if (!SetupInner_bool) {

        Nelem_inner = Nelem_;
        FourierOrder_inner = FourierOrder_;
        ElemOrder_inner = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated;
        if (Xc_inner.Dim() == Nelem_ * ElemOrder_ * 3) {
            std::cout << "Cursory check on Xc_inner size, doesn't guarantee nodes are distributed in panel-based Cheb quadrature!" << std::endl;
            Xc_updated = Xc_inner;
            r_updated = r_inner;
        } else {
            std::cout << "Xc_inner does not match expected number of nodes using panel-based quadrature; updating Xc_inner." << std::endl;
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            // TODO: NEED TO CHANGE IN THE FUTURE, BUT FOR NOW: assume r = r(x), no y and z dependency.
            // collect src_x and trg_x
            sctl::Vector<Real> src_x(Xc_inner.Dim()/3);
            sctl::Vector<Real> trg_x(Xc_updated.Dim()/3);
            for (sctl::Long node_ind = 0; node_ind < Xc_inner.Dim()/3; node_ind++) {
                src_x[node_ind] = Xc_inner[node_ind*3];
            }
            for (sctl::Long node_ind = 0; node_ind < Xc_updated.Dim()/3; node_ind++) {
                trg_x[node_ind] = Xc_updated[node_ind*3];
            }
            InterpR(r_updated, r_inner, src_x, trg_x);
        }
        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_inner_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated);
        elem_lst_inner = elem_lst_inner_;

        SetupInner_bool = true;
        Xc_inner = Xc_updated;
        r_inner = r_updated;
        min_radius = sctl::omp_par::reduce(r_inner.begin(), r_inner.Dim());
        elem_lst_inner.GetNodeCoord(&X_inner, &Xn_inner, nullptr);
    } else {
        // This clause is only for when called from user function; shouldn't be here from Setup().
        std::cout << "Note: Inner already setup, nothing done." << std::endl;
    }
    return std::make_tuple(Xc_inner, r_inner);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> Annular<Real>::SetupOuter(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_) {
    if (!SetupOuter_bool) {
        Nelem_outer = Nelem_;
        FourierOrder_outer = FourierOrder_;
        ElemOrder_outer = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated;
        if (Xc_outer.Dim() == Nelem_ * ElemOrder_ * 3) {
            std::cout << "Cursory check on Xc_outer size, doesn't guarantee nodes are distributed in panel-based Cheb quadrature!" << std::endl;
            Xc_updated = Xc_outer;
            r_updated = r_outer;
            // TODO: will add adaptive quadrature later. Also check what kind of adaptive quadrature is already implemented in CSBQ, should be some Fourier zooming in..?
        } else {
            std::cout << "Xc_outer does not match expected number of nodes using panel-based quadrature; updating Xc_inner." << std::endl;
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            // TODO: NEED TO CHANGE IN THE FUTURE, BUT FOR NOW: assume r = r(x), no y and z dependency.
            // collect src_x and trg_x
            sctl::Vector<Real> src_x(Xc_inner.Dim()/3);
            sctl::Vector<Real> trg_x(Xc_updated.Dim()/3);
            for (sctl::Long node_ind = 0; node_ind < Xc_inner.Dim()/3; node_ind++) {
                src_x[node_ind] = Xc_inner[node_ind*3];
            }
            for (sctl::Long node_ind = 0; node_ind < Xc_updated.Dim()/3; node_ind++) {
                trg_x[node_ind] = Xc_updated[node_ind*3];
            }
            InterpR(r_updated, r_inner, src_x, trg_x);
        }
        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_outer_(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated); // SlenderElemList will compute outward normal, so must be use Xn_outer to get correct orientation.
        elem_lst_outer = elem_lst_outer_;

        SetupOuter_bool = true;
        Xc_outer = Xc_updated;
        r_outer = r_updated;
        elem_lst_outer.GetNodeCoord(&X_outer, &Xn_outer, nullptr);
        // Correct normal of outer channel to point into fluid domain.
        Xn_outer *= -1.;
    } else {
        std::cout << "Note: Outer already setup, nothing done." << std::endl;
    }
    return std::make_tuple(Xc_outer, r_outer);
}

template <class Real>
void Annular<Real>::GetInnerCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out) {
    // DBC: if not setup, causes assertion error.
    SCTL_ASSERT(SetupInner_bool); 
    (*X_out) = X_inner;
    if (Xn_out) {
        (*Xn_out) = Xn_inner;
    }
}

template <class Real>
void Annular<Real>::GetOuterCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out) {
    SCTL_ASSERT(SetupOuter_bool); 
    (*X_out) = X_outer;
    if (Xn_out) {
        (*Xn_out) = Xn_outer;
    }
}

template <class Real>
void Annular<Real>::GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn) {
    SCTL_ASSERT(SetupInner_bool && SetupOuter_bool); 
    if (Xn) {
        X->ReInit(X_inner.Dim() + X_outer.Dim());
        Xn->ReInit(Xn_inner.Dim() + Xn_outer.Dim());
        for (sctl::Long ind=0; ind<X_inner.Dim(); ind++) {
            (*X)[ind] = X_inner[ind];
            (*Xn)[ind] = Xn_inner[ind];
        }
        sctl::Long offset = X_inner.Dim();
        for (sctl::Long ind=0; ind<X_outer.Dim(); ind++) {
            (*X)[ind+offset] = X_outer[ind];
            (*Xn)[ind+offset] = Xn_outer[ind];
        }
    } else {
        // no need to return Xn
        X->ReInit(X_inner.Dim() + X_outer.Dim());
        for (sctl::Long ind=0; ind<X_inner.Dim(); ind++) {
            (*X)[ind] = X_inner[ind];
        }
        sctl::Long offset = X_inner.Dim();
        for (sctl::Long ind=0; ind<X_outer.Dim(); ind++) {
            (*X)[ind+offset] = X_outer[ind];
        }
    }
    
}

// x in [0,1], reinitializes x only, set y and z to 0.5... 
template <class Real>
sctl::Vector<Real> Annular<Real>::GetCenterLine(sctl::Long Nelem_, sctl::Long ElemOrder_) {
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
Real Annular<Real>::GetMinRadius() {
    SCTL_ASSERT(SetupInner_bool);
    return min_radius;
}

template <class Real>
sctl::Vector<sctl::Long> Annular<Real>::InDomain(sctl::Vector<Real> X_) {
    sctl::Long Nnodes = X_.Dim()/3;
    sctl::Vector<sctl::Long> in_domain(Nnodes);
    for (sctl::Long nd = 0; nd < Nnodes; nd++) {
        const Real x = X_[nd*3+0];
        const Real y = X_[nd*3+1];
        const Real z = X_[nd*3+2];
        if (x < Xc_inner[0] || x < Xc_outer[0] || x > Xc_inner[Xc_inner.Dim()-1] || x > Xc_outer[Xc_outer.Dim()-1]) {
            // x value already out of channel, don't consider.
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
            };
            // find panel x belongs to
            sctl::Long x_inner_ind = search_panel(x, Nelem_inner, ElemOrder_inner, Xc_inner);
            sctl::Long x_outer_ind = search_panel(x, Nelem_outer, ElemOrder_outer, Xc_outer);
            // interpolate r to this x value
            sctl::Vector<Real> r_inner_interp(1);
            sctl::Vector<Real> src_r(ElemOrder_inner, (sctl::Iterator<Real>)r_inner.begin()+x_inner_ind*ElemOrder_inner, false);
            sctl::Vector<Real> src_x(ElemOrder_inner*3, (sctl::Iterator<Real>)Xc_inner.begin() + x_inner_ind*ElemOrder_inner*3, false);
            sctl::Vector<Real> trg_x(3, (sctl::Iterator<Real>)X_.begin() + 3*nd, false);
            InterpR(r_inner_interp, src_r, src_x, trg_x);
            
            sctl::Vector<Real> r_outer_interp(1);
            sctl::Vector<Real> src_r2(ElemOrder_outer, (sctl::Iterator<Real>)r_outer.begin()+x_outer_ind*ElemOrder_outer, false);
            sctl::Vector<Real> src_x2(ElemOrder_outer*3, (sctl::Iterator<Real>)Xc_outer.begin() + x_outer_ind*ElemOrder_outer*3, false);
            InterpR(r_outer_interp, src_r2, src_x2, trg_x);
            // compare r^2 with (y^2+z^2)
            in_domain[nd] = ((r_outer_interp[0]*r_outer_interp[0] - y*y - z*z) > 1e-5) && (y*y + z*z - (r_inner_interp[0]*r_inner_interp[0]) > 1e-5);
        }
    }
    return in_domain;
}

template class Annular<float>;
template class Annular<double>;
