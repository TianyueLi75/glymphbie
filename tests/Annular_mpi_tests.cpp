#include <GlymphBIE.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

// TESTS FOR MPI construction only -- will be comparing with serial version, so assumes serial code works fine.
// TODO: (later) when Nelem not divisible by Nmpi

template <class Real>
void test_centerline(sctl::Comm comm)
{
    std::cout << "\n Starting test_centerline test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    // Cheb from serial code
    auto Xc_serial = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    // Cheb from mpi
    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);

    // std::cout << "DEBUGGING, Check size. Xc mpi size on rank " << comm.Rank() << " is " << Xc_mpi.Dim() << ", constructor says elem count is " << cnt << ", displacement is " << dsp << "; size of serial is " << Xc_serial.Dim() << std::endl;

    SCTL_ASSERT(Xc_mpi.Dim() == 3 * cnt * ElemOrder);

    sctl::Long offset = dsp*ElemOrder*3;
    for (sctl::Long ind=0; ind < Xc_mpi.Dim(); ind++) {
        // std::cout << "DEBUGGING, Rank " << comm.Rank() << ", mpi at ind " << ind << " is " << Xc_mpi[ind] << ", serial looks at ind " << ind+offset << ", has value " << Xc_serial[offset+ind] << std::endl;
        SCTL_ASSERT(sctl::fabs(Xc_mpi[ind] - Xc_serial[offset + ind]) < 1e-6);
    }

    std::cout << "Passed test_centerline test. \n" << std::endl;
}

template <class Real>
void straight_getnodes(sctl::Comm comm) 
{
    std::cout << "\n Starting get_nodes test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    auto Xc_serial = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(Xc_serial.Dim()/3);
    r1 = 0.3;
    sctl::Vector<Real> r2(Xc_serial.Dim()/3);
    r2 = 0.5;
    sctl::Vector<Real> drdx(Xc_serial.Dim()/3); // no outward velocity on wall.
    drdx = 0.;
    Annular<Real> straight(Xc_serial,Xc_serial,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx); 

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    // std::cout << "on rank " << comm.Rank() << ", cnt is " << cnt << ", dsp is " << dsp << ", size of Xc_mpi here is " << Xc_mpi.Dim() << std::endl;
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdx_mpi(Xc_mpi.Dim()/3); // no outward velocity on wall.
    drdx_mpi = 0.;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdx_mpi,drdx_mpi);
    
    Real rad_out = straight.GetMinRadius();
    // MPI Reduce to collect min of all min radii on each process
    sctl::Vector<Real> rad_loc(1);
    rad_loc[0] = rad_out;
    sctl::Vector<Real> rad_all(1);
    rad_all[0] = 0.;
    comm.Allreduce((sctl::Iterator<Real>) rad_loc.begin(), (sctl::Iterator<Real>) rad_all.begin(), 1, sctl::CommOp::MIN);
    // Assert on each process that the collected min radii equals the expected min radius.
    SCTL_ASSERT(sctl::fabs(rad_all[0] - r1[0]) < 1e-6);
    
    // Check coordinates of channel setup using MPI match serial setup.
    sctl::Vector<Real> X_inner, X_outer, X_all, Xn_inner, Xn_outer, Xn_all;
    straight.GetInnerCoord(&X_inner, &Xn_inner);
    straight.GetOuterCoord(&X_outer, &Xn_outer);
    straight.GetNodeCoord(&X_all, &Xn_all);

    sctl::Vector<Real> X_inner_mpi, X_outer_mpi, X_all_mpi, Xn_inner_mpi, Xn_outer_mpi, Xn_all_mpi;
    straight_mpi.GetInnerCoord(&X_inner_mpi, &Xn_inner_mpi);
    straight_mpi.GetOuterCoord(&X_outer_mpi, &Xn_outer_mpi);
    straight_mpi.GetNodeCoord(&X_all_mpi, &Xn_all_mpi);

    // std::cout << "rank " << comm.Rank() << ", overall Nelem is " << Nelem << ", ElemOrder = " << ElemOrder << ", FO = " << FourierOrder <<  "; size of Xinner mpi is " << X_inner_mpi.Dim() << ",  Xinner total is " << X_inner.Dim() << std::endl;

    SCTL_ASSERT(X_inner_mpi.Dim() == 3*cnt*ElemOrder*FourierOrder);

    sctl::Long offset = dsp*ElemOrder*3*FourierOrder;
    for (sctl::Long ind=0; ind < X_inner_mpi.Dim(); ind++) {
        SCTL_ASSERT(sctl::fabs(X_inner_mpi[ind] - X_inner[offset + ind]) < 1e-6);
    }
    for (sctl::Long ind=0; ind < X_outer_mpi.Dim(); ind++) {
        SCTL_ASSERT(sctl::fabs(X_outer_mpi[ind] - X_outer[offset + ind]) < 1e-6);
    }

    std::cout << "Passed get_nodes test. \n" << std::endl;
}

template <class Real>
void in_domain(sctl::Comm comm) 
{
    std::cout << "\n Starting in_domain test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdx_mpi(Xc_mpi.Dim()/3); // no outward velocity on wall.
    drdx_mpi = 0.;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder, drdx_mpi, drdx_mpi);

    sctl::Vector<Real> Xtrg(3);
    Xtrg = {0.75,0.75,0.75};
    sctl::Vector<sctl::Long> in_dom = straight_mpi.InDomain(Xtrg);
    std::cout << "On rank " << comm.Rank() << "; in domain is " << in_dom[0] << std::endl;

    // Collect in_domain boolean from all process
    sctl::Vector<sctl::Long> in_dom_loc = in_dom;
    sctl::Vector<sctl::Long> in_dom_all(in_dom.Dim());
    in_dom_all = 0;
    comm.Allreduce((sctl::Iterator<Real>) in_dom_loc.begin(), (sctl::Iterator<Real>) in_dom_all.begin(), in_dom.Dim(), sctl::CommOp::MAX); // will be 1 if ==1 on any process; =0 if all ==0.
    in_dom = in_dom_all;

    SCTL_ASSERT(in_dom.Dim()==1);
    SCTL_ASSERT(in_dom[0]);

    std::cout << "Passed in_domain test. \n" << std::endl;
}

template <class Real>
void get_vslip(sctl::Comm comm) 
{
    std::cout << "\n Starting get_vslip test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;

    auto Xc_serial = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(Xc_serial.Dim()/3);
    r1 = 0.3;
    sctl::Vector<Real> r2(Xc_serial.Dim()/3);
    r2 = 0.5;
    sctl::Vector<Real> drdx(Xc_serial.Dim()/3); // no outward velocity on wall.
    drdx = 0.5;
    Annular<Real> straight(Xc_serial,Xc_serial,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx); 
    sctl::Vector<Real> v2 = straight.GetVslip();

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    // std::cout << "on rank " << comm.Rank() << ", cnt is " << cnt << ", dsp is " << dsp << ", size of Xc_mpi here is " << Xc_mpi.Dim() << std::endl;
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdx_mpi(Xc_mpi.Dim()/3); 
    drdx_mpi = 0.5;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder, drdx_mpi, drdx_mpi);
    sctl::Vector<Real> v2_mpi = straight_mpi.GetVslip_mpi();

    sctl::Long offset = dsp*ElemOrder*3*FourierOrder;
    for (sctl::Long ind=0; ind < v2_mpi.Dim(); ind++) {
        SCTL_ASSERT(sctl::fabs(v2_mpi[ind] - v2[offset + ind]) < 1e-6);
    }

    std::cout << "Passed get_vslip test. \n" << std::endl;
}



int main(int argc, char** argv)
{
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real=double;

    {
        sctl::Comm comm = sctl::Comm::World();
        test_centerline<Real>(comm);
        straight_getnodes<Real>(comm);
        in_domain<Real>(comm);
        get_vslip<Real>(comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}