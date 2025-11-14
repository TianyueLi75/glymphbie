#ifndef _ANNULAR_HPP_
#define _ANNULAR_HPP_

#include <csbq.hpp>

template <class Real> class Annular {

    public:
        // TODO: add MPI support for communicator later.
        Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_);
        ~Annular();

        // TODO: keeping these functions public for now for user to update geometry easily.
        void SetInnerXc(const sctl::Vector<Real> X_in);
        void SetOuterXc(const sctl::Vector<Real> X_in);
        void SetInnerR(const sctl::Vector<Real> r_in);
        void SetOuterR(const sctl::Vector<Real> r_in);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Setup(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> SetupInner(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_); // TODO: Compute aspect ratio in inner setup since this will be the most limiting factor
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>> SetupOuter(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_);
        void GetInnerCoord(sctl::Vector<Real>* X_out);
        void GetOuterCoord(sctl::Vector<Real>* X_out);
        void GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn); // TODO: enable comm support later.
        static sctl::Vector<Real> GetCenterLine(const sctl::Long Nelem_, const sctl::Long ElemOrder_);
        Real GetMinRadius();

        // TODO: add read/write and plotting support

        // TODO: setup adaptive grid, and return new centerline after each solve.



    private:
        sctl::Vector<Real> Xc_inner, Xc_outer, r_inner, r_outer; // input values
        sctl::Vector<Real> X_inner, X_outer, Xn_inner, Xn_outer; // computed values from CSBQ
        Real min_radius; // computed value
        sctl::Long Nelem_inner, Nelem_outer, ElemOrder_inner, ElemOrder_outer, FourierOrder_inner, FourierOrder_outer; 

        sctl::SlenderElemList<Real> elem_lst_outer, elem_lst_inner;

        bool SetupInner_bool = false;
        bool SetupOuter_bool = false;

        bool CheckCenterLine(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const bool check_inner);
        void InterpR(sctl::Vector<Real>& trg_r, const sctl::Vector<Real> src_r, const sctl::Vector<Real> src_x, const sctl::Vector<Real> trg_x);
};










#endif