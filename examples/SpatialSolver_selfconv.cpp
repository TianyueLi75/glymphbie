// Verification for the spatial solver -- 
// Self convergence.
#include <GlymphBIE.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept> 

// Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
template <class Real> 
sctl::Vector<Real> bg_poiseuille(const Real dpdx, const Real mu, const sctl::Vector<Real> X)
{
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    U.SetZero();
    for (sctl::Long i = 0; i < N; i++) {
        sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
        U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4./mu; 
    }
    return U;
}

// For a sinusoidal channel on the outside and circular channel on the inside, 
// First perform self-convergence tests at random points in domain
template <class Real> 
void concentric_sine_selfconv(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm, const sctl::Long write_ref)
{

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    const Real R_in = 0.13;
    const Real R_out = 0.48;

    // // Channel 1: inner channel wavy, outer fixed
    // const auto getRin = [&R_in](const Real x) {
    //     return R_in * (0.3+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    // };
    // const auto getRout = [&R_out](const Real x) {
    //     return R_out;
    // };

    // Channel 2: inner channel fixed, outer wavy
    const auto getRin = [&R_in](const Real x) {
        return R_in;
    };
    const auto getRout = [&R_out](const Real x) {
        return R_out * (0.7+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    };

    // Get outlet radius info for flux calculations.
    Real R_in_x1 = getRin(1.);
    Real R_out_x1 = getRout(1.);

    // Make Annular channel
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    for (sctl::Long i=0; i<r1.Dim(); i++) {
        r1[i] = getRin(Xc[i*3+0]);
    }
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    for (sctl::Long i=0; i<r2.Dim(); i++) {
        r2[i] = getRout(Xc[i*3+0]);
    }
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> sinusoid(Xc,Xc,r1,r2);
    sinusoid.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);

    // Make StokesBIO
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    Real S_scal = 1.0;
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    LPO.AddElemList(sinusoid.GetInnerElemList(), "inner");
    LPO.AddElemList(sinusoid.GetOuterElemList(), "outer");
    LPO.SetAccuracy(tol);
    sctl::Vector<Real> X_annular;
    sinusoid.GetNodeCoord(&X_annular, nullptr);
    LPO.SetTargetCoord(X_annular);
    LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
    sctl::Vector<Real> NormalOrient(X_annular.Dim());
    NormalOrient = 1.;
    sctl::Long size_inner = Nelem*ElemOrder*FourierOrder*3;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }
    // Periodic single-layer nullspace: far-field quadrature weights and total
    // surface area, used to remove the surface-mean density during the solve.
    sctl::SlenderElemList<Real> elem_inner = sinusoid.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = sinusoid.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

    sctl::Vector<Real> X1temp, X2temp;
    sinusoid.GetInnerCoord(&X1temp, nullptr);
    sinusoid.GetOuterCoord(&X2temp, nullptr);
    std::cout << "size of X1 temp: " << X1temp.Dim() << ", size of X2 temp: " << X2temp.Dim() << std::endl;
    sinusoid.WriteVTK("../vis/sine_channel_inner", "../vis/sine_channel_outer", X1temp, X2temp, comm);

    // Create lambda function
    const auto BIO = [&LPO, D_scal, &NormalOrient, &elem_inner, &elem_outer, size_inner, &wts_inner, &wts_outer, surface_area, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
        // Remove the surface-mean density (periodic single-layer nullspace).
        sctl::Vector<Real> sigma_mean = ComputeSigmaMean<Real>(sigma, elem_inner, elem_outer, size_inner, wts_inner, wts_outer, surface_area, comm);
        sctl::Vector<Real> sigma0 = sigma;
        AddConstVec<Real>(sigma0, sigma_mean*(Real)(-1));
        U->SetZero();
        LPO.ComputePotential(*U, sigma0);
        if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma0*D_scal*NormalOrient; // double-layer jump condition on surface.
        AddConstVec<Real>(*U, sigma_mean); // add back surface-mean density.
    };

    // Setup gmres
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    // Utot = Ubg + Uwall = 0 on wall. 
    sctl::Vector<Real> vslip = sinusoid.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in_x1, &R_out_x1, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out_x1 - R_in_x1) + R_in_x1;
        sctl::Vector<Real> r_wts = wts * (R_out_x1 - R_in_x1);

        sctl::Vector<Real> theta_nodes(theta_ord);
        const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;
        for (sctl::Long i=0; i<theta_ord; i++) {
            theta_nodes[i] = theta_wt * i;
        }

        // form mesh grid of x-y points, same-r-all-thetas order.
        sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Xtrg[3*node_ind + 0] = 1.0; // evaluation at x=1.0 only
                Xtrg[3*node_ind + 1] = 0.5 + r_nodes[i] * sctl::cos<Real>(theta_nodes[j]); 
                Xtrg[3*node_ind + 2] = 0.5 + r_nodes[i] * sctl::sin<Real>(theta_nodes[j]); 
            }
        }

        // Compute U at meshgrid:
        sctl::Vector<Real> Utrg;
        LPO.SetTargetCoord(Xtrg);
        BIO(&Utrg, sigma);
        Utrg += bg_poiseuille(dpdx,mu,Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        // Real Qsa = 0.;
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * r_nodes[i] * theta_wt * r_wts[i];
            }
        }
        Real Q2 = Q * mu / (-dpdx) / R_in_x1 / R_in_x1 / R_in_x1 / R_in_x1; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(45,64);
    std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder <<  "; Flux Q computed = " << Q << std::endl; // should be the same for all processes

    // Set target points in between channels 
    sctl::Long Nelem_trg = 1;
    sctl::Long ElemOrder_trg = 10;
    sctl::Long FourierOrder_trg = 16;
    sctl::Vector<Real> Xc_trg = Annular<Real>::GetCenterLine(Nelem_trg, ElemOrder_trg);
    sctl::Vector<Real> r1_trg(ElemOrder_trg*Nelem_trg);
    for (sctl::Long i=0; i<r1_trg.Dim(); i++) {
        r1_trg[i] = getRin(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> r2_trg(ElemOrder_trg*Nelem_trg);
    for (sctl::Long i=0; i<r2_trg.Dim(); i++) {
        r2_trg[i] = getRout(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> drdt_trg(ElemOrder_trg*Nelem_trg); // no outward velocity on wall.
    drdt_trg = 0.;
    sctl::Vector<Real> r3 = r1_trg + 1./3. * (r2_trg - r1_trg);
    sctl::Vector<Real> r4 = r1_trg + 2./3. * (r2_trg - r1_trg);
    Annular<Real> channel_trg(Xc_trg,Xc_trg,r3,r4);
    channel_trg.Setup(Nelem_trg, ElemOrder_trg, FourierOrder_trg,drdt_trg,drdt_trg);
    sctl::Vector<Real> X_inner, X_outer;
    channel_trg.GetInnerCoord(&X_inner, nullptr);
    channel_trg.GetOuterCoord(&X_outer, nullptr);
    sctl::Vector<Real> U_inner, U_outer;
    LPO.SetTargetCoord(X_inner);
    BIO(&U_inner, sigma);
    LPO.SetTargetCoord(X_outer);
    BIO(&U_outer, sigma);
    U_inner += bg_poiseuille(dpdx,mu,X_inner);
    U_outer += bg_poiseuille(dpdx,mu,X_outer);

    // write to file or read and compare for error. 
    std::string filename_in = "../out/SelfConv/Concentric_sin_4_32_U_exact_"+std::to_string(comm.Rank())+"_inner.txt";
    std::string filename_out = "../out/SelfConv/Concentric_sin_4_32_U_exact_"+std::to_string(comm.Rank())+"_outer.txt";
    std::string filename_q = "../out/SelfConv/Concentric_sin_4_32_Q_exact_"+std::to_string(comm.Rank())+".txt";
    if (write_ref==1) {
        U_inner.Write(filename_in.c_str());
        U_outer.Write(filename_out.c_str());
        std::ofstream outfile(filename_q.c_str());
        outfile << Q << std::endl;
        outfile.close();
        // visualization
        channel_trg.WriteVTK("../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_inner","../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_outer", U_inner, U_outer, comm);
    } else {
        sctl::Vector<Real> U_ref_in, U_ref_out;
        U_ref_in.Read(filename_in.c_str());
        U_ref_out.Read(filename_out.c_str());
        Real Q_ref;
        std::ifstream infile(filename_q.c_str());
        infile >> Q_ref;
        infile.close();

        // Combine error and max_u values on inner and outer channels
        sctl::Vector<Real> err = U_inner - U_ref_in;
        Real max_err = 0;
        for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        Real max_u = 0.;
        for (const auto e : U_ref_in) max_u = std::max<Real>(max_u, sctl::fabs(e));
        err = U_outer - U_ref_out;
        for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        for (const auto e : U_ref_out) max_u = std::max<Real>(max_u, sctl::fabs(e));

        Real Qerr = fabs(Q_ref - Q) / fabs(Q_ref);

        sctl::Vector<Real> err_loc(1);
        err_loc[0] = max_err;
        sctl::Vector<Real> err_all(1);
        err_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);
        
        sctl::Vector<Real> u_loc(1);
        u_loc[0] = max_u;
        sctl::Vector<Real> u_all(1);
        u_all[0] = 0.;
        comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

        if (!comm.Rank()) {
            std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; max error = " << std::setprecision(10) << err_all[0] << "; max relative error = " << err_all[0] / u_all[0] << "; relative error of flux Q is " << Qerr << std::endl;
        }
    }

}

template <class Real> 
void concentric_sine_selfconv_mpi(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm, const sctl::Long write_ref)
{

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    const Real R_in = 0.13;
    const Real R_out = 0.48;

    // // Channel 1: inner channel wavy, outer fixed
    // const auto getRin = [&R_in](const Real x) {
    //     return R_in * (0.3+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    // };
    // const auto getRout = [&R_out](const Real x) {
    //     return R_out;
    // };

    // Channel 2: inner channel fixed, outer wavy
    const auto getRin = [&R_in](const Real x) {
        return R_in;
    };
    const auto getRout = [&R_out](const Real x) {
        return R_out * (0.7+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    };

    // Get outlet radius info for flux calculations.
    Real R_in_x1 = getRin(1.);
    Real R_out_x1 = getRout(1.);

    // Make Annular channel
    auto [cnt, dsp, Xc] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    sctl::Vector<Real> r1(Xc.Dim()/3);
    for (sctl::Long i=0; i<r1.Dim(); i++) {
        r1[i] = getRin(Xc[i*3+0]);
    }
    sctl::Vector<Real> r2(Xc.Dim()/3);
    for (sctl::Long i=0; i<r2.Dim(); i++) {
        r2[i] = getRout(Xc[i*3+0]);
    }
    sctl::Vector<Real> drdt(Xc.Dim()/3); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> sinusoid(Xc,Xc,r1,r2, comm);
    sinusoid.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdt,drdt);

    // Make StokesBIO
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    Real S_scal = 1.0;
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    LPO.AddElemList(sinusoid.GetInnerElemList(), "inner");
    LPO.AddElemList(sinusoid.GetOuterElemList(), "outer");
    LPO.SetAccuracy(tol);
    sctl::Vector<Real> X_annular;
    sinusoid.GetNodeCoord(&X_annular, nullptr);
    LPO.SetTargetCoord(X_annular);
    LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
    sctl::Vector<Real> NormalOrient(X_annular.Dim());
    NormalOrient = 1.;
    sctl::Long size_inner = Xc.Dim() * FourierOrder;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }
    // Periodic single-layer nullspace: far-field quadrature weights and total
    // surface area, used to remove the surface-mean density during the solve.
    sctl::SlenderElemList<Real> elem_inner = sinusoid.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = sinusoid.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

    sctl::Vector<Real> X1temp, X2temp;
    sinusoid.GetInnerCoord(&X1temp, nullptr);
    sinusoid.GetOuterCoord(&X2temp, nullptr);
    sinusoid.WriteVTK("../vis/sine_channel_inner", "../vis/sine_channel_outer", X1temp, X2temp, comm);

    // Create lambda function
    const auto BIO = [&LPO, D_scal, &NormalOrient, &elem_inner, &elem_outer, size_inner, &wts_inner, &wts_outer, surface_area, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
        // Remove the surface-mean density (periodic single-layer nullspace).
        sctl::Vector<Real> sigma_mean = ComputeSigmaMean<Real>(sigma, elem_inner, elem_outer, size_inner, wts_inner, wts_outer, surface_area, comm);
        sctl::Vector<Real> sigma0 = sigma;
        AddConstVec<Real>(sigma0, sigma_mean*(Real)(-1));
        U->SetZero();
        LPO.ComputePotential(*U, sigma0);
        if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma0*D_scal*NormalOrient; // double-layer jump condition on surface.
        AddConstVec<Real>(*U, sigma_mean); // add back surface-mean density.
    };

    // Setup gmres
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    // Utot = Ubg + Uwall = 0 on wall. 
    sctl::Vector<Real> vslip = sinusoid.GetVslip_mpi(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in_x1, &R_out_x1, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out_x1 - R_in_x1) + R_in_x1;
        sctl::Vector<Real> r_wts = wts * (R_out_x1 - R_in_x1);

        sctl::Vector<Real> theta_nodes(theta_ord);
        const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;
        for (sctl::Long i=0; i<theta_ord; i++) {
            theta_nodes[i] = theta_wt * i;
        }

        // form mesh grid of x-y points, same-r-all-thetas order.
        sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Xtrg[3*node_ind + 0] = 1.0; // evaluation at x=1.0 only
                Xtrg[3*node_ind + 1] = 0.5 + r_nodes[i] * sctl::cos<Real>(theta_nodes[j]); 
                Xtrg[3*node_ind + 2] = 0.5 + r_nodes[i] * sctl::sin<Real>(theta_nodes[j]); 
            }
        }

        // Compute U at meshgrid:
        sctl::Vector<Real> Utrg;
        LPO.SetTargetCoord(Xtrg);
        BIO(&Utrg, sigma);
        Utrg += bg_poiseuille(dpdx,mu,Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        // Real Qsa = 0.;
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * r_nodes[i] * theta_wt * r_wts[i];
                // Qsa += r_nodes[i] * theta_wt * r_wts[i];
            }
        }
        
        Real Q2 = Q * mu / (-dpdx) / R_in_x1 / R_in_x1 / R_in_x1 / R_in_x1; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(45,64);
    if (!comm.Rank()) {
        std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder <<  "; Flux Q computed = " << Q << std::endl; // should be the same for all processes
    }

    // Set target points in between channels 
    // For now: same targets on all process.
    sctl::Long Nelem_trg = 1;
    sctl::Long ElemOrder_trg = 10;
    sctl::Long FourierOrder_trg = 16;
    sctl::Vector<Real> Xc_trg = Annular<Real>::GetCenterLine(Nelem_trg, ElemOrder_trg);
    sctl::Vector<Real> r1_trg(ElemOrder_trg*Nelem_trg);
    for (sctl::Long i=0; i<r1_trg.Dim(); i++) {
        r1_trg[i] = getRin(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> r2_trg(ElemOrder_trg*Nelem_trg);
    for (sctl::Long i=0; i<r2_trg.Dim(); i++) {
        r2_trg[i] = getRout(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> drdt_trg(ElemOrder_trg*Nelem_trg); // no outward velocity on wall.
    drdt_trg = 0.;
    sctl::Vector<Real> r3 = r1_trg + 1./3. * (r2_trg - r1_trg);
    sctl::Vector<Real> r4 = r1_trg + 2./3. * (r2_trg - r1_trg);
    Annular<Real> channel_trg(Xc_trg,Xc_trg,r3,r4);
    channel_trg.Setup(Nelem_trg, ElemOrder_trg, FourierOrder_trg,drdt_trg,drdt_trg);
    sctl::Vector<Real> X_inner, X_outer;
    channel_trg.GetInnerCoord(&X_inner, nullptr);
    channel_trg.GetOuterCoord(&X_outer, nullptr);
    sctl::Vector<Real> U_inner, U_outer;
    LPO.SetTargetCoord(X_inner);
    BIO(&U_inner, sigma);
    LPO.SetTargetCoord(X_outer);
    BIO(&U_outer, sigma);
    U_inner += bg_poiseuille(dpdx,mu,X_inner);
    U_outer += bg_poiseuille(dpdx,mu,X_outer);

    // write to file or read and compare for error. 
    std::string filename_in = "../out/SelfConv/Concentric_sin_4_32_U_exact_0_inner.txt";
    std::string filename_out = "../out/SelfConv/Concentric_sin_4_32_U_exact_0_outer.txt";
    std::string filename_q = "../out/SelfConv/Concentric_sin_4_32_Q_exact_0.txt";
    if (write_ref==1) {
        if (!comm.Rank()) {
            U_inner.Write(filename_in.c_str());
            U_outer.Write(filename_out.c_str());
            std::ofstream outfile(filename_q.c_str());
            outfile << Q << std::endl;
            outfile.close();
            // visualization
            channel_trg.WriteVTK("../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_inner","../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_outer", U_inner, U_outer, comm.Self());
        }
    } else {
        if (!comm.Rank()) {
            sctl::Vector<Real> U_ref_in, U_ref_out;
            U_ref_in.Read(filename_in.c_str());
            U_ref_out.Read(filename_out.c_str());
            Real Q_ref;
            std::ifstream infile(filename_q.c_str());
            infile >> Q_ref;
            infile.close();

            // Combine error and max_u values on inner and outer channels
            sctl::Vector<Real> err = U_inner - U_ref_in;
            Real max_err = 0;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            Real max_u = 0.;
            for (const auto e : U_ref_in) max_u = std::max<Real>(max_u, sctl::fabs(e));
            err = U_outer - U_ref_out;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            for (const auto e : U_ref_out) max_u = std::max<Real>(max_u, sctl::fabs(e));

            Real Qerr = fabs(Q_ref - Q) / fabs(Q_ref);

            std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; max error = " << std::setprecision(10) << max_err << "; max relative error = " << max_err / max_u << "; relative error of flux Q is " << Qerr << std::endl;

        }
    }

}



int main(int argc, char** argv)
{
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real=double;

    {
        sctl::Comm comm = sctl::Comm::World();
        // Read parameters from input
        long Nelem = std::stol(argv[1]); // number of elements
        long ElemOrder = std::stol(argv[2]); // Cheb order per element
        long FourierOrder = std::stol(argv[3]);  // number of Fourier nodes
        // Set default physical and problem parameters
        double mu = 1.0;
        double dpdx = -1.;

        if (Nelem==4 && FourierOrder == 32) {
            concentric_sine_selfconv<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 1);
            // concentric_sine_selfconv_mpi<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 1);
        } else {
            concentric_sine_selfconv<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 0);
            // concentric_sine_selfconv_mpi<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 0);
        }
        
        
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

