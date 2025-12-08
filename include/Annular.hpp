#ifndef _ANNULAR_HPP_
#define _ANNULAR_HPP_

// #include "csbq.hpp"
#include "csbq/slender_element.hpp"
#include "csbq/slender_element.cpp"

template <class Real> class Annular {

    public:
        Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_);
        Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_, sctl::Comm comm_);
        ~Annular();

        // TODO: keeping these functions public for now for user to update geometry easily.
        void SetInnerXc(const sctl::Vector<Real> X_in);
        void SetOuterXc(const sctl::Vector<Real> X_in);
        void SetInnerR(const sctl::Vector<Real> r_in);
        void SetOuterR(const sctl::Vector<Real> r_in);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Setup(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Setup_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, sctl::Vector<Real> drdt_outer_, bool force_setup = false);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupInner(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, bool force_setup = false);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupOuter(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupInner_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, bool force_setup = false);
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupOuter_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        void GetInnerCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out);
        void GetOuterCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out);
        void GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn); // TODO: enable comm support later.
        static sctl::Vector<Real> GetCenterLine(const sctl::Long Nelem_, const sctl::Long ElemOrder_);
        static std::tuple<sctl::Long, sctl::Long, sctl::Vector<Real>> GetCenterLine_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, sctl::Comm comm);
        Real GetMinRadius();
        sctl::Vector<sctl::Long> InDomain(sctl::Vector<Real> X_); // given a vector of x positions, return which nodes are in bewteen outer and inner channels V.S. not.
        sctl::Vector<sctl::Long> InDomain_mpi(sctl::Vector<Real> X_);
        sctl::SlenderElemList<Real> GetInnerElemList();
        sctl::SlenderElemList<Real> GetOuterElemList();
        sctl::Vector<Real> GetVslip();
        sctl::Vector<Real> GetVslip_mpi();
        sctl::Vector<Real> GetVslipInner();
        sctl::Vector<Real> GetVslipOuter();
        sctl::Vector<Real> GetVslipInner_mpi();
        sctl::Vector<Real> GetVslipOuter_mpi();

        void WriteVTK(std::string filename_inner, std::string filename_outer, sctl::Vector<Real> F_inner, sctl::Vector<Real> F_outer, sctl::Comm comm);

        // TODO: setup adaptive grid, and return new centerline after each solve.



    private:
        sctl::Vector<Real> Xc_inner, Xc_outer, r_inner, r_outer; // input values
        sctl::Vector<Real> X_inner, X_outer, Xn_inner, Xn_outer; // computed values from CSBQ; NOTE: Normal points into fluid.
        Real min_radius = 0.; // computed value
        sctl::Long Nelem_inner, Nelem_outer, ElemOrder_inner, ElemOrder_outer, FourierOrder_inner, FourierOrder_outer; 
        sctl::Vector<Real> drdt_inner, drdt_outer;

        // MPI specific params
        sctl::Comm comm;
        sctl::Long loc_elem_dsp_inner, loc_elem_dsp_outer;
        sctl::Long loc_elem_cnt_inner, loc_elem_cnt_outer;

        sctl::SlenderElemList<Real> elem_lst_outer, elem_lst_inner;

        bool SetupInner_bool = false;
        bool SetupOuter_bool = false;

        void InterpR(sctl::Vector<Real>& trg_r, const sctl::Vector<Real> src_r, const sctl::Vector<Real> src_x, const sctl::Vector<Real> trg_x);
};









#endif