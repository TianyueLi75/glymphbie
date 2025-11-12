#include "Annular.hpp"

// TODO: input Xc may be z-periodic instead of x-, should add a conversion function here.
template <class Real>
Annular<Real>::Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_) {
    SCTL_ASSERT(Xc_inner_.Dim()==r_inner_.Dim());
    SCTL_ASSERT(Xc_outer_.Dim()==r_outer_.Dim());
    sctl::Vector<Real> Xc_inner = Xc_inner_;
    sctl::Vector<Real> Xc_outer = Xc_outer_;
    sctl::Vector<Real> r_inner = r_inner_;
    sctl::Vector<Real> r_outer = r_outer_;
    sctl::Vector<Real> X_inner, X_outer, Xn_inner, Xn_outer;
    Real aspect_ratio = 0.;
}

template <class Real>
Annular<Real>::~Annular() {} // TODO: This is it for the descrutor?


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
        auto [Xc_outer_updated, r_outer_updated] = SetupOuter(Nelem_, ElemOrder_, FourierOrder_);
        Xc_inner_updated = Xc_inner;
        r_inner_updated = r_inner;
    } else if (SetupOuter_bool) {
        // Outer set up, inner not
        auto [Xc_inner_updated, r_inner_updated] = SetupInner(Nelem_, ElemOrder_, FourierOrder_);
        Xc_outer_updated = Xc_outer;
        r_outer_updated = r_outer;
    } else {
        // Neither is set up
        auto [Xc_inner_updated, r_inner_updated] = SetupInner(Nelem_, ElemOrder_, FourierOrder_);
        auto [Xc_outer_updated, r_outer_updated] = SetupOuter(Nelem_, ElemOrder_, FourierOrder_);
    }    
    return std::make_tuple(Xc_inner_updated, Xc_outer_updated, r_inner_updated, r_outer_updated);
}

//TODO: may need to add functionality for scanning to add a vector of ElemOrder, etc. Right now assumes all panels have the same orders in GL and fourier.
template <class Real>
bool Annular<Real>::CheckCenterLine(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const bool check_inner) {
    sctl::Long expected_len = Nelem_ * ELemOrder_;
    if (check_inner) {
        return Xc_inner.Dim() == expected_len;
    } else {
        return Xc_outer.Dim() == expected_len;
    }
}


template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> Annular<Real>::SetupInner(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_) {
    if (!SetupInner_bool) {
        sctl::Long Nelem_inner = Nelem_;
        sctl::Long FourierOrder_inner = FourierOrder_;
        sctl::Long ElemOrder_inner = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated;
        if (CheckCenterLine(Nelem_, ElemOrder_, true)) {
            std::cout << "Cursory check on Xc_inner size, doesn't guarantee nodes are distributed in panel-based GL quadrature!" << std::endl;
            Xc_updated = Xc_inner;
            r_updated = r_inner;
        } else {
            std::cout << "Xc_inner does not match expected number of nodes using panel-based quadrature; updating Xc_inner." << std::endl;
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            // TODO: make interpolation matrix for r

            // TODO: get r_updated;
            // r_updated = ; 
        }
        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_inner(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated);

        SetupInner_bool = true;
        Xc_inner = Xc_updated;
        r_inner = r_updated;
        // TODO: compute aspect ratio
    } else {
        // This clause is only for when called from user function; shouldn't be here from Setup().
        std::cout << "Note: Inner already setup, nothing done." << std::endl;
    }
    return std::make_tuple(Xc_inner, r_inner);
}

template <class Real>
std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> Annular<Real>::SetupOuter(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_) {
    if (!SetupOuter_bool) {
        sctl::Long Nelem_outer = Nelem_;
        sctl::Long FourierOrder_outer = FourierOrder_;
        sctl::Long ElemOrder_outer = ElemOrder_;

        sctl::Vector<Real> Xc_updated, r_updated;
        if (CheckCenterLine(Nelem_, ElemOrder_, true)) {
            std::cout << "Cursory check on Xc_inner size, doesn't guarantee nodes are distributed in panel-based GL quadrature!" << std::endl;
            Xc_updated = Xc_outer;
            r_updated = r_outer;
            // TODO: will add adaptive quadrature later. Also check what kind of adaptive quadrature is already implemented in CSBQ, should be some Fourier zooming in..?
        } else {
            std::cout << "Xc_inner does not match expected number of nodes using panel-based quadrature; updating Xc_inner." << std::endl;
            Xc_updated = GetCenterLine(Nelem_, ElemOrder_);
            // TODO: make interpolation matrix for r

            // TODO: get r_updated;
            // r_updated = ; 
        }
        // TODO: add MPI functions.
        sctl::Vector<sctl::Long> ElemOrderVec(Nelem_);
        sctl::Vector<sctl::Long> FourierOrderVec(Nelem_);
        ElemOrderVec = ElemOrder_;
        FourierOrderVec = FourierOrder_;
        sctl::SlenderElemList<Real> elem_lst_outer(ElemOrderVec, FourierOrderVec, Xc_updated, r_updated); // TODO: check whether without orient, surface normal is computed correctly.

        SetupInner_bool = true;
        Xc_outer = Xc_updated;
        r_outer = r_updated;
    } else {
        std::cout << "Note: Outer already setup, nothing done." << std::endl;
    }
    return std::make_tuple(Xc_outer, r_outer);
}

template <class Real>
void Annular<Real>::GetInnerCoord(sctl::Vector<Real>* X_out) {
    // DBC: if not setup, causes assertion error.
    SCTL_ASSERT(SetupInner_bool); 
    X_out = X_inner;
}

template <class Real>
void Annular<Real>::GetOuterCoord(sctl::Vector<Real>* X_out) {
    SCTL_ASSERT(SetupOuter_bool); 
    X_out = X_outer;
}

template <class Real>
void Annular<Real>::GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn) {
    SCTL_ASSERT(SetupInner_bool && SetupOuter_bool); 
    if (Xn) {
        X.ReInit(X_inner.Dim() + X_outer.Dim());
        Xn.ReInit(Xn_inner.Dim() + Xn_outer.Dim());
        for (sctl::Long ind=0; ind<X_inner.Dim(); ind++) {
            X[ind] = X_inner[ind];
            Xn[ind] = Xn_inner[ind];
        }
        sctl::Long offset = X_inner.Dim();
        for (sctl::Long ind=0; ind<X_outer.Dim(); ind++) {
            X[ind+offset] = X_outer[ind];
            Xn[ind+offset] = Xn_outer[ind];
        }
    } else {
        // no need to return Xn
        X.ReInit(X_inner.Dim() + X_outer.Dim());
        for (sctl::Long ind=0; ind<X_inner.Dim(); ind++) {
            X[ind] = X_inner[ind];
        }
        sctl::Long offset = X_inner.Dim();
        for (sctl::Long ind=0; ind<X_outer.Dim(); ind++) {
            X[ind+offset] = X_outer[ind];
        }
    }
    
}

// x in [0,1], reinitializes x only, set y and z to 0.5... 
// TODO: may need to interpolate x and y? Can we start with just straight channel?
template <class Real>
sctl::Vector<Real> Annular<Real>::GetCenterLine(sctl::Long Nelem_, sctl::Long ElemOrder_) {
    // TODO: returns centerline nodes when x=[0,1] divided into <Nelem> panels of <ElemOrder> number of GL nodes each
    const sctl::Vector<Real> nodes_01 = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder_);
    SCTL_ASSERT(ElemOrder_ == nodes_01.Dim());
    sctl::Vector<Real> Xc(ElemOrder_ * Nelem_);
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
Real Annular<Real>::GetAspectRatio() {
    SCTL_ASSERT(SetupInner_bool);
    return aspect_ratio;
}

template class Annular<float>;
template class Annular<double>;
