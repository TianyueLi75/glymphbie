// Verification for the spatial solver -- Compare with exact solutions, manufactured solutions, and potentially FDM or other published code.
#include <GlymphBIE.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

template <class Real> 
void concentric_poiseuille(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
{
    const Real R_in = 0.3;
    const Real R_out = 0.48;

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    // Make Annular channel
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    r1 = R_in;
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    r2 = R_out;
    sctl::Vector<Real> drdx(ElemOrder*Nelem); // no outward velocity on wall.
    drdx = 0.;
    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx);
    // std::cout << "DEBUG: Xc: "<< std::endl;
    // for (int i=0; i<Xc.Dim(); i++) {
    //     std::cout << Xc[i] << std::endl;
    // }

    // Function to get exact flow at location X given center Xc and other params.
    // NOTE: This assumes Xvec taken as an intermediate channel with the same parameters for simplicity. Otherwise an interpolation is needed.
    const auto uexact = [R_in, R_out, dpdx, mu, &Xc, Nelem, ElemOrder, FourierOrder](const sctl::Vector<Real> Xvec) {
        sctl::Long N = Xvec.Dim()/3;
        SCTL_ASSERT(Xvec.Dim() == Nelem * ElemOrder * FourierOrder * 3);
        sctl::Vector<Real> Uvec(Xvec.Dim());
        Uvec.SetZero();
        // compute constants first -- TODO: arithemetic order to reduce numerical error
        Real coeff = 1./4./mu;
        // Get rmax
        Real rmax_num = R_out*R_out - R_in*R_in;
        Real Rratio = R_out / R_in;
        Real rmax_denom = 2*sctl::log<Real>(Rratio);
        Real rmax2 = rmax_num / rmax_denom;
        // Compute u = 1/4/mu*dpdx*(R_out^2-r^2-2rmax^2ln*(R_out/r))
        for (sctl::Long i=0; i<Nelem; i++) {
            for (sctl::Long j=0; j<ElemOrder; j++) {
                // get Xc at this index
                sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>)Xc.begin() + (i*ElemOrder+j)*3, false);
                for (sctl::Long k=0; k<FourierOrder; k++) {
                    sctl::Long starting_ind = (i*ElemOrder*FourierOrder + j*FourierOrder + k)*3;
                    sctl::Vector<Real> Xhere(3, (sctl::Iterator<Real>) Xvec.begin() + starting_ind, false);
                    sctl::Vector<Real> Rhere = Xhere - Xc_here;
                    Real r2 = Rhere[1]*Rhere[1] + Rhere[2]*Rhere[2];
                    Real rratio = R_out / sctl::sqrt<Real>(r2);
                    Real u = -1.*coeff * dpdx * (R_out*R_out - r2 - 2. * rmax2 * sctl::log<Real>(rratio));
                    Uvec[starting_ind + 0] = u; // flow always just in the x direction.
                }
            }
        }
        return Uvec;
    };

    // Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
    const auto bg_poiseuille = [dpdx, mu](const sctl::Vector<Real> X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        U.SetZero();
        for (sctl::Long i = 0; i < N; i++) {
            sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
            U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4./mu; 
        }
        return U;
    };

    // Make StokesBIO
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    Real S_scal = 1.0;
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    LPO.AddElemList(straight.GetInnerElemList(), "inner");
    LPO.AddElemList(straight.GetOuterElemList(), "outer");
    LPO.SetAccuracy(tol);
    sctl::Vector<Real> X_annular;
    straight.GetNodeCoord(&X_annular, nullptr);
    LPO.SetTargetCoord(X_annular);
    LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
    sctl::Vector<Real> NormalOrient(X_annular.Dim());
    NormalOrient = 1.;
    sctl::Long size_inner = Nelem*ElemOrder*FourierOrder*3;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    // std::cout << "size inner = " << size_inner << ", size outer = " << size_outer << std::endl;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }

    // Create lambda function
    const auto BIO = [&LPO, D_scal, &NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
        U->SetZero();
        LPO.ComputePotential(*U, sigma);
        if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma*D_scal*NormalOrient;
    };

    // Setup gmres
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    // Utot = Ubg + Uwall = 0 on wall. 
    sctl::Vector<Real> vslip = straight.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(X_annular);
    // std::cout << "DEBUG, bg poiseuille values, should be equal in magnitude for all nodes on each channel." << std::endl;
    for (int i=0; i<X_annular.Dim()/3; i++) {
        Real magv2 = vbg[i*3+0]*vbg[i*3+0] + vbg[i*3+1]*vbg[i*3+1] + vbg[i*3+2]*vbg[i*3+2];
        // std::cout << "vector = " << vbg[i*3+0] << ", " << vbg[i*3+1] << ", " << vbg[i*3+2] << ", magnitude = " << magv2 << std::endl;
        // std::cout << "vslip (shoudl be 0): " << vslip[i*3+0] << ", " << vslip[i*3+1] << ", " << vslip[i*3+2] << std::endl;
    }

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compare to exact, report error.
    // Create two cylindrical channels of different radius to test U at different r. (NOTE: Can do sinusoidal channel later)
    sctl::Vector<Real> r3(ElemOrder*Nelem);
    r3 = R_in + 1./3. * (R_out - R_in);
    sctl::Vector<Real> r4(ElemOrder*Nelem);
    r4 = R_in + 2./3. * (R_out - R_in);

    Annular<Real> channel_trg(Xc,Xc,r3,r4);
    channel_trg.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx);
    sctl::Vector<Real> X_inner, X_outer;
    channel_trg.GetInnerCoord(&X_inner, nullptr);
    channel_trg.GetOuterCoord(&X_outer, nullptr);
    sctl::Vector<Real> U_inner, U_outer;
    LPO.SetTargetCoord(X_inner);
    BIO(&U_inner, sigma);
    LPO.SetTargetCoord(X_outer);
    BIO(&U_outer, sigma);
    U_inner += bg_poiseuille(X_inner);
    U_outer += bg_poiseuille(X_outer);
    sctl::Vector<Real> Uexact_inner = uexact(X_inner);
    sctl::Vector<Real> Uexact_outer = uexact(X_outer);
    Real max_err_inner = 0.;
    Real max_err_outer = 0.;
    sctl::Vector<Real> Diff_inner = U_inner - Uexact_inner;
    sctl::Vector<Real> Diff_outer = U_outer - Uexact_outer;
    // std::cout << " ======================= DEBUG U eval." << std::endl;
    // Real mag_exact2 = Uexact_inner[0]*Uexact_inner[0] + Uexact_inner[1]*Uexact_inner[1] + Uexact_inner[2]*Uexact_inner[2];
    // std::cout << "uexact inner = " << Uexact_inner[0] << ", " << Uexact_inner[1] << ", " << Uexact_inner[2] << ", mag = " << mag_exact2 << std::endl;
    // Real mag_solve2 = U_inner[0]*U_inner[0] + U_inner[1]*U_inner[1] + U_inner[2]*U_inner[2];
    // std::cout << "usolve inner = " << U_inner[0] << ", " << U_inner[1] << ", " << U_inner[2] << ", mag = " << mag_solve2 << std::endl;
    // mag_exact2 = Uexact_outer[0]*Uexact_outer[0] + Uexact_outer[1]*Uexact_outer[1] + Uexact_outer[2]*Uexact_outer[2];
    // std::cout << "uexact outer = " << Uexact_outer[0] << ", " << Uexact_outer[1] << ", " << Uexact_outer[2] << ", mag = " << mag_exact2 << std::endl;
    // mag_solve2 = U_outer[0]*U_outer[0] + U_outer[1]*U_outer[1] + U_outer[2]*U_outer[2];
    // std::cout << "usolve outer = " << U_outer[0] << ", " << U_outer[1] << ", " << U_outer[2] << ", mag = " << mag_solve2 << std::endl;
    for (auto e : Diff_inner) max_err_inner = std::max<Real>(max_err_inner, sctl::fabs(e));
    for (auto e : Diff_outer) max_err_outer = std::max<Real>(max_err_outer, sctl::fabs(e));
    Real avg_max_err = 0.5*(max_err_inner + max_err_outer);
    std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Averaged max error at two different r: "<< avg_max_err << std::endl;

    // Check flux calculation against analytical solutions.
    const Real alpha = R_out / R_in;
    const Real beta = alpha;
    const Real K = alpha*beta - 1.0;
    const Real Q_ana = sctl::const_pi<Real>()/8. *((K+1.)*(K+1.) - 1.0 - 2.*K*K/sctl::log<Real>(K+1.));

    // // Check parameters to Tithof paper analytical U
    // Real rho1 = (R_in + 1./3. * (R_out - R_in)) / R_in;
    // Real rho2 = (R_in + 2./3. * (R_out - R_in)) / R_in;
    // std::cout << rho1 << std::endl;
    // Real U_T1 = 1/4. * (alpha*alpha - rho1*rho1 - (alpha*alpha-1)*sctl::log<Real>(alpha/rho1)/sctl::log<Real>(alpha));
    // Real U_T2 = 1/4. * (alpha*alpha - rho2*rho2 - (alpha*alpha-1)*sctl::log<Real>(alpha/rho2)/sctl::log<Real>(alpha));
    // std::cout << "debug U value, U eval here is " << U_T1 << ", U paper scaled by r_1^2: " << U_T1 * R_in * R_in * (-1 * dpdx) << std::endl;
    // U_T1 *= R_in * R_in * (-1 * dpdx) / mu;
    // U_T2 *= R_in * R_in * (-1 * dpdx) / mu;
    
    // std::cout << "Difference between U in paper and U here are " << fabs(U_T1 - Uexact_inner[0]) << ", " << fabs(U_T2 - Uexact_outer[0]) << std::endl;

    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in, &R_out, dpdx, mu, &sigma, &BIO, &bg_poiseuille](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out - R_in) + R_in;
        sctl::Vector<Real> r_wts = wts * (R_out - R_in);

        sctl::Vector<Real> theta_nodes(theta_ord);
        const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;
        for (int i=0; i<theta_ord; i++) {
            theta_nodes[i] = theta_wt * i;
        }

        // form mesh grid of x-y points, same-r-all-thetas order.
        sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
        for (int i=0; i<r_ord; i++) {
            for (int j=0; j<theta_ord; j++) {
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
        Utrg += bg_poiseuille(Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        for (int i=0; i<r_ord; i++) {
            for (int j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * r_nodes[i] * theta_wt * r_wts[i];
            }
        }

        Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(20,64);

    std::cout << "Flux difference: Q analytical = " << std::setprecision(10) << Q_ana << ", Q computed = " << Q << "; rel error is " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;

}

// template <class Real> 
// void concentric_poiseuille_mpi(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
// {
//     const Real R_in = 0.3;
//     const Real R_out = 0.48;

//     Real gmres_tol = 1e-11;
//     Real tol = 1e-14;
//     if (FourierOrder < 16) {
//         gmres_tol = 1e-9;
//     }

//     // Make Annular channel
//     auto [cnt, dsp, Xc] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
//     sctl::Vector<Real> r1(Xc.Dim()/3);
//     r1 = R_in;
//     sctl::Vector<Real> r2(Xc.Dim()/3);
//     r2 = R_out;
//     sctl::Vector<Real> drdx(Xc.Dim()/3); // no outward velocity on wall.
//     drdx = 0.;
//     Annular<Real> straight(Xc,Xc,r1,r2, comm);
//     straight.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdx,drdx);
//     // std::cout << "DEBUG: Xc: "<< std::endl;
//     // for (int i=0; i<Xc.Dim(); i++) {
//     //     std::cout << Xc[i] << std::endl;
//     // }

//     // Function to get exact flow at location X given center Xc and other params.
//     // NOTE: This assumes Xvec taken as an intermediate channel with the same parameters for simplicity. Otherwise an interpolation is needed.
//     const auto uexact = [R_in, R_out, dpdx, mu, &Xc, Nelem, ElemOrder, FourierOrder](const sctl::Vector<Real> Xvec) {
//         sctl::Long N = Xvec.Dim()/3;
//         SCTL_ASSERT(Xvec.Dim() == Nelem * ElemOrder * FourierOrder * 3);
//         sctl::Vector<Real> Uvec(Xvec.Dim());
//         Uvec.SetZero();
//         // compute constants first -- TODO: arithemetic order to reduce numerical error
//         Real coeff = 1./4./mu;
//         // Get rmax
//         Real rmax_num = R_out*R_out - R_in*R_in;
//         Real Rratio = R_out / R_in;
//         Real rmax_denom = 2*sctl::log<Real>(Rratio);
//         Real rmax2 = rmax_num / rmax_denom;
//         // Compute u = 1/4/mu*dpdx*(R_out^2-r^2-2rmax^2ln*(R_out/r))
//         for (sctl::Long i=0; i<Nelem; i++) {
//             for (sctl::Long j=0; j<ElemOrder; j++) {
//                 // get Xc at this index
//                 sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>)Xc.begin() + (i*ElemOrder+j)*3, false);
//                 for (sctl::Long k=0; k<FourierOrder; k++) {
//                     sctl::Long starting_ind = (i*ElemOrder*FourierOrder + j*FourierOrder + k)*3;
//                     sctl::Vector<Real> Xhere(3, (sctl::Iterator<Real>) Xvec.begin() + starting_ind, false);
//                     sctl::Vector<Real> Rhere = Xhere - Xc_here;
//                     Real r2 = Rhere[1]*Rhere[1] + Rhere[2]*Rhere[2];
//                     Real rratio = R_out / sctl::sqrt<Real>(r2);
//                     Real u = -1.*coeff * dpdx * (R_out*R_out - r2 - 2. * rmax2 * sctl::log<Real>(rratio));
//                     Uvec[starting_ind + 0] = u; // flow always just in the x direction.
//                 }
//             }
//         }
//         return Uvec;
//     };

//     // Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
//     const auto bg_poiseuille = [dpdx, mu](const sctl::Vector<Real> X) {
//         const sctl::Long N = X.Dim()/3;
//         sctl::Vector<Real> U(N*3);
//         U.SetZero();
//         for (sctl::Long i = 0; i < N; i++) {
//             sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
//             U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4./mu; 
//         }
//         return U;
//     };

//     // Make StokesBIO
//     // Real min_rad = straight.GetMinRadius();
//     // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
//     Real S_scal = 1.0;
//     Real D_scal = 1.0;
//     StokesBIO<Real> LPO(S_scal, D_scal, comm);
//     LPO.AddElemList(straight.GetInnerElemList(), "inner");
//     LPO.AddElemList(straight.GetOuterElemList(), "outer");
//     LPO.SetAccuracy(tol);
//     sctl::Vector<Real> X_annular;
//     straight.GetNodeCoord(&X_annular, nullptr);
//     LPO.SetTargetCoord(X_annular);
//     LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
//     sctl::Vector<Real> NormalOrient(X_annular.Dim());
//     NormalOrient = 1.;
//     sctl::Long size_inner = Nelem*ElemOrder*FourierOrder*3;
//     sctl::Long size_outer = X_annular.Dim() - size_inner;
//     // std::cout << "size inner = " << size_inner << ", size outer = " << size_outer << std::endl;
//     for (sctl::Long ind=0; ind<size_outer; ind++) {
//         // Normal of outer channel element list points outward by default.
//         NormalOrient[size_inner + ind] = -1.;
//     }

//     // Create lambda function
//     const auto BIO = [&LPO, D_scal, &NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
//         U->SetZero();
//         LPO.ComputePotential(*U, sigma);
//         if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma*D_scal*NormalOrient;
//     };

//     // Setup gmres
//     sctl::GMRES<Real> solver(comm);
//     sctl::KrylovPrecond<Real> ksp;
//     sctl::Vector<Real> sigma;
//     // Utot = Ubg + Uwall = 0 on wall. 
//     sctl::Vector<Real> vslip = straight.GetVslip_mpi(); // no EXTRA wall velocity
//     sctl::Vector<Real> vbg = bg_poiseuille(X_annular);
//     // std::cout << "DEBUG, bg poiseuille values, should be equal in magnitude for all nodes on each channel." << std::endl;
//     for (int i=0; i<X_annular.Dim()/3; i++) {
//         Real magv2 = vbg[i*3+0]*vbg[i*3+0] + vbg[i*3+1]*vbg[i*3+1] + vbg[i*3+2]*vbg[i*3+2];
//         // std::cout << "vector = " << vbg[i*3+0] << ", " << vbg[i*3+1] << ", " << vbg[i*3+2] << ", magnitude = " << magv2 << std::endl;
//         // std::cout << "vslip (shoudl be 0): " << vslip[i*3+0] << ", " << vslip[i*3+1] << ", " << vslip[i*3+2] << std::endl;
//     }

//     // Solve gmres
//     solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

//     // Compare to exact, report error.
//     // Create two cylindrical channels of different radius to test U at different r. 
//     sctl::Vector<Real> r3(Xc.Dim()/3);
//     r3 = R_in + 1./3. * (R_out - R_in);
//     sctl::Vector<Real> r4(Xc.Dim()/3);
//     r4 = R_in + 2./3. * (R_out - R_in);

//     Annular<Real> channel_trg(Xc,Xc,r3,r4,comm);
//     channel_trg.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdx,drdx);
//     sctl::Vector<Real> X_inner, X_outer;
//     channel_trg.GetInnerCoord(&X_inner, nullptr);
//     channel_trg.GetOuterCoord(&X_outer, nullptr);
//     sctl::Vector<Real> U_inner, U_outer;
//     LPO.SetTargetCoord(X_inner);
//     BIO(&U_inner, sigma);
//     LPO.SetTargetCoord(X_outer);
//     BIO(&U_outer, sigma);
//     U_inner += bg_poiseuille(X_inner);
//     U_outer += bg_poiseuille(X_outer);
//     sctl::Vector<Real> Uexact_inner = uexact(X_inner);
//     sctl::Vector<Real> Uexact_outer = uexact(X_outer);
//     Real max_err_inner = 0.;
//     Real max_err_outer = 0.;
//     sctl::Vector<Real> Diff_inner = U_inner - Uexact_inner;
//     sctl::Vector<Real> Diff_outer = U_outer - Uexact_outer;
//     for (auto e : Diff_inner) max_err_inner = std::max<Real>(max_err_inner, sctl::fabs(e));
//     for (auto e : Diff_outer) max_err_outer = std::max<Real>(max_err_outer, sctl::fabs(e));
//     Real avg_max_err = 0.5*(max_err_inner + max_err_outer);
//     // MPI Reduce to collect min of all min radii on each process
//     sctl::Vector<Real> err_loc(1);
//     err_loc[0] = avg_max_err;
//     sctl::Vector<Real> err_all(1);
//     err_all[0] = 0.;
//     comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);

//     std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Averaged max error at two different r: "<< err_all[0]<< std::endl;

//     // Check flux calculation against analytical solutions.
//     const Real alpha = R_out / R_in;
//     const Real beta = alpha;
//     const Real K = alpha*beta - 1.0;
//     const Real Q_ana = sctl::const_pi<Real>()/8. *((K+1.)*(K+1.) - 1.0 - 2.*K*K/sctl::log<Real>(K+1.));

//     // // Check parameters to Tithof paper analytical U
//     // Real rho1 = (R_in + 1./3. * (R_out - R_in)) / R_in;
//     // Real rho2 = (R_in + 2./3. * (R_out - R_in)) / R_in;
//     // std::cout << rho1 << std::endl;
//     // Real U_T1 = 1/4. * (alpha*alpha - rho1*rho1 - (alpha*alpha-1)*sctl::log<Real>(alpha/rho1)/sctl::log<Real>(alpha));
//     // Real U_T2 = 1/4. * (alpha*alpha - rho2*rho2 - (alpha*alpha-1)*sctl::log<Real>(alpha/rho2)/sctl::log<Real>(alpha));
//     // std::cout << "debug U value, U eval here is " << U_T1 << ", U paper scaled by r_1^2: " << U_T1 * R_in * R_in * (-1 * dpdx) << std::endl;
//     // U_T1 *= R_in * R_in * (-1 * dpdx) / mu;
//     // U_T2 *= R_in * R_in * (-1 * dpdx) / mu;
    
//     // std::cout << "Difference between U in paper and U here are " << fabs(U_T1 - Uexact_inner[0]) << ", " << fabs(U_T2 - Uexact_outer[0]) << std::endl;

//     // Compute Q by doing integral over outlet flow
//     const auto getQ = [&LPO, &R_in, &R_out, dpdx, mu, &sigma, &BIO, &bg_poiseuille](const sctl::Long r_ord, const sctl::Long theta_ord) {
//         SCTL_ASSERT(dpdx<0.);
        
//         sctl::Vector<Real> nds, wts;
//         sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
//         sctl::Vector<Real> r_nodes = nds * (R_out - R_in) + R_in;
//         sctl::Vector<Real> r_wts = wts * (R_out - R_in);

//         sctl::Vector<Real> theta_nodes(theta_ord);
//         const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;
//         for (int i=0; i<theta_ord; i++) {
//             theta_nodes[i] = theta_wt * i;
//         }

//         // form mesh grid of x-y points, same-r-all-thetas order.
//         sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
//         for (int i=0; i<r_ord; i++) {
//             for (int j=0; j<theta_ord; j++) {
//                 sctl::Long node_ind = i*theta_ord + j;
//                 Xtrg[3*node_ind + 0] = 1.0; // evaluation at x=1.0 only
//                 Xtrg[3*node_ind + 1] = 0.5 + r_nodes[i] * sctl::cos<Real>(theta_nodes[j]); 
//                 Xtrg[3*node_ind + 2] = 0.5 + r_nodes[i] * sctl::sin<Real>(theta_nodes[j]); 
//             }
//         }

//         // Compute U at meshgrid:
//         sctl::Vector<Real> Utrg;
//         LPO.SetTargetCoord(Xtrg);
//         BIO(&Utrg, sigma);
//         Utrg += bg_poiseuille(Xtrg);

//         // Quadrature to integrate
//         Real Q = 0.;
//         for (int i=0; i<r_ord; i++) {
//             for (int j=0; j<theta_ord; j++) {
//                 sctl::Long node_ind = i*theta_ord + j;
//                 Real u = Utrg[3*node_ind + 0];
//                 Q += u * r_nodes[i] * theta_wt * r_wts[i];
//             }
//         }

//         Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA

//         return Q2;

//     };

//     Real Q = getQ(20,64); // For now just compute flux on all processes, can also check whether all return same value.

//     std::cout << "Flux difference: Q analytical = " << std::setprecision(10) << Q_ana << ", Q computed = " << Q << "; rel error is " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;

// }

// /*
// For eccentric channels, use fixed centerlines and radii to prevent intersection.
template <class Real> 
void eccentric_poiseuille(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
{

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    const Real R_in = 0.13;
    const Real R_out = 0.48;
    const Real Xc_y_out = 0.5;
    const Real Xc_z_out = 0.5;
    const Real Xc_y_in = 0.6;
    const Real Xc_z_in = 0.4;

    // Compute parameters for analytical solutions
    const Real alpha = R_out / R_in;
    const Real beta = alpha;
    const Real K = alpha*beta - 1.0;
    const Real c = sctl::sqrt<Real>((Xc_y_out-Xc_y_in)*(Xc_y_out-Xc_y_in) + (Xc_z_out-Xc_z_in)*(Xc_z_out-Xc_z_in));
    // std::cout << "check c, should be 0.1414: " << c << std::endl;

    const Real alpha2 = alpha * alpha;
    const Real alpha4 = alpha2 * alpha2;
    const Real eccen = c / R_in;
    const Real eccen2 = eccen * eccen;
    const Real F = (alpha2 - 1 + eccen2)/2./eccen;
    const Real M = sctl::sqrt<Real>((F*F-alpha2));
    const Real A = 0.5*sctl::log<Real>((F+M)/(F-M));
    const Real B = 0.5*sctl::log<Real>((F-eccen+M)/(F-eccen-M));
    // std::cout << "alpha (should be 3.69230769231) = " << alpha << ", K (should be 12.6331360947) = " << K << ", F (should be 6.35036282335) = " << F << ", A (should be 1.13757114107) = " << A << ", B (should be 2.3446025812) = " << B << std::endl;
    const Real inf_sum = 0.0205399624783; // converged finite sum at N~20 in Desmos, specific to these parameters.
    const Real Q_ana = sctl::const_pi<Real>()/8. *((alpha4-1.) - 4*eccen2*M*M/(B-A) - 8 * eccen2 * M * M * inf_sum);

    // Make Annular channel
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> Xc_inner = Xc;
    for (sctl::Long i=0; i<Xc_inner.Dim()/3; i++) {
        Xc_inner[i*3+1] = Xc_y_in;
        Xc_inner[i*3+2] = Xc_z_in;
    }
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    r1 = R_in;
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    r2 = R_out;
    sctl::Vector<Real> drdx(ElemOrder*Nelem); // no outward velocity on wall.
    drdx = 0.;
    Annular<Real> straight(Xc_inner,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx);
    // std::cout << "DEBUG: Xc: "<< std::endl;
    // for (int i=0; i<Xc.Dim(); i++) {
    //     std::cout << Xc[i] << std::endl;
    // }

    // Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
    const auto bg_poiseuille = [dpdx, mu](const sctl::Vector<Real> X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        U.SetZero();
        for (sctl::Long i = 0; i < N; i++) {
            sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
            U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4. / mu; // TODO: check direction
        }
        return U;
    };

    // Make StokesBIO
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    Real S_scal = 1.0;
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    LPO.AddElemList(straight.GetInnerElemList(), "inner");
    LPO.AddElemList(straight.GetOuterElemList(), "outer");
    LPO.SetAccuracy(tol);
    sctl::Vector<Real> X_annular;
    straight.GetNodeCoord(&X_annular, nullptr);
    LPO.SetTargetCoord(X_annular);
    LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
    sctl::Vector<Real> NormalOrient(X_annular.Dim());
    NormalOrient = 1.;
    sctl::Long size_inner = Nelem*ElemOrder*FourierOrder*3;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    // std::cout << "size inner = " << size_inner << ", size outer = " << size_outer << std::endl;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }

    // Create lambda function
    const auto BIO = [&LPO, D_scal, &NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
        U->SetZero();
        LPO.ComputePotential(*U, sigma);
        if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma*D_scal*NormalOrient;
    };

    // Setup gmres
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    // Utot = Ubg + Uwall = 0 on wall. 
    sctl::Vector<Real> vslip = straight.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    // r1 and r2 differ at different theta's!
    const auto getQ = [&LPO, &R_in, &R_out, &Xc_y_in, &Xc_y_out, &Xc_z_in, &Xc_z_out, dpdx, mu, &sigma, &BIO, &bg_poiseuille](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]

        const auto positiveRoot = [](double a, double b, double c) {
            double disc = b*b - 4*a*c;
            if(disc < 0) {
                std::cout << "discreminant in r_max finding is negative, fails." << std::endl;
                return 0.0; // no real solution
            }

            // std::cout << "b = " << b << ", c = " << c << ", disc = " << disc << ", sqrt(disc) = " << sctl::sqrt<Real>(disc) << ", std sqrt: " << std::sqrt(disc) << std::endl;
            
            double r1 = (-b + std::sqrt(disc)) / (2.*a);
            double r2 = (-b - std::sqrt(disc)) / (2.*a);
            double rmax = std::max(r1,r2);
            // std::cout << "max of two roots is " << rmax << ", min of two roots is " << r1+r2-rmax << std::endl;
            return rmax;
        };

        // form mesh grid of x-y points, same-theta-all-r order.
        sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
        sctl::Vector<Real> r_wts_trg(r_ord * theta_ord);
        sctl::Vector<Real> r_nds_trg(r_ord * theta_ord);
        sctl::Vector<Real> r_out_trg(theta_ord);
        const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;

        Real Q_eq = 0.;
        Real d1 = Xc_y_in-Xc_y_out;
        Real d2 = Xc_z_in-Xc_z_out;
        Real d2sq = d1*d1+d2*d2;

        for (int i=0; i<theta_ord; i++) {
            Real theta_val = theta_wt * i;
            // Real dy = R_out*sctl::cos<Real>(theta_val) + Xc_y_out - Xc_y_in;
            // Real dz = R_out*sctl::sin<Real>(theta_val) + Xc_z_out - Xc_z_in;
            // Real r_out = sctl::sqrt<Real>((dz*dz) + (dy*dy));

            Real cosT = sctl::cos<Real>(theta_val);
            Real sinT = sctl::sin<Real>(theta_val);
            double b = 2.*(d1*cosT + d2*sinT);
            // std::cout << "theta = " << theta_val << ", cosT = " << cosT << ", sinT = " << sinT << std::endl;
            double c = d2sq - R_out*R_out;
            double rho_max = positiveRoot(1.0, b, c);

            Real r_out = rho_max; // DEBUGGING: whether the Cheb integral didn't work purely becuase of r_out.

            sctl::Vector<Real> r_nodes = nds * (r_out - R_in) + R_in;
            sctl::Vector<Real> r_wts = wts * (r_out - R_in);
            
            // Real drho = (rho_max - R_in) / r_ord;

            r_out_trg[i] = r_out;

            for (int j=0; j<r_ord; j++) {
                sctl::Long node_ind = i*r_ord + j;
                Xtrg[3*node_ind + 0] = 1.0; // evaluation at x=1.0 only
                Xtrg[3*node_ind + 1] = Xc_y_in + r_nodes[j] * sctl::cos<Real>(theta_val); 
                Xtrg[3*node_ind + 2] = Xc_z_in + r_nodes[j] * sctl::sin<Real>(theta_val); 

                r_wts_trg[node_ind] = r_wts[j];
                r_nds_trg[node_ind] = r_nodes[j];

                // // DEBUG equidist
                // double rho = R_in + (j + 0.5) * drho; // midpoint
                // double x = Xc_y_in + rho * cosT;
                // double y = Xc_z_in + rho * sinT;
                // Q_eq += rho * drho * theta_wt;
                // Xtrg[3*node_ind + 1] = x;
                // Xtrg[3*node_ind + 2] = y;
                // // std::cout << "X,Y from equidist: " << x <<", " << y <<std::endl;
                // r_wts_trg[node_ind] = drho;
                // r_nds_trg[node_ind] = rho;
            }
        }
        // std::cout << "Q surface area from equidist points is " << Q_eq << "; SA is " << sctl::const_pi<Real>() * (R_out*R_out - R_in*R_in) << std::endl;

        // std::cout << "Debugging surface integral." << std::endl;
        // Real Q1 = 0.;
        // for (int i=0; i<theta_ord; i++) {
        //     for (int j=0; j<r_ord; j++) {
        //         sctl::Long node_ind = i*r_ord + j;
        //         Q1 += theta_wt * r_wts_trg[node_ind] * r_nds_trg[node_ind];
        //     }
        // }
        // std::cout << "Q1 from quadrature is " << Q1 << "; SA is " << sctl::const_pi<Real>() * (R_out*R_out - R_in*R_in) << std::endl;

        // Compute U at meshgrid:
        sctl::Vector<Real> Utrg;
        LPO.SetTargetCoord(Xtrg);
        BIO(&Utrg, sigma);
        Utrg += bg_poiseuille(Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        for (int i=0; i<theta_ord; i++) {
            for (int j=0; j<r_ord; j++) {
                sctl::Long node_ind = i*r_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * theta_wt * r_wts_trg[node_ind] * r_nds_trg[node_ind];
            }
        }

        Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(45,64);
    std::cout <<"Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Flux relative error = " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;
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
    for (int i=0; i<r1.Dim(); i++) {
        r1[i] = getRin(Xc[i*3+0]);
    }
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    for (int i=0; i<r2.Dim(); i++) {
        r2[i] = getRout(Xc[i*3+0]);
    }
    sctl::Vector<Real> drdx(ElemOrder*Nelem); // no outward velocity on wall.
    drdx = 0.;
    Annular<Real> sinusoid(Xc,Xc,r1,r2);
    sinusoid.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx);

    // Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
    const auto bg_poiseuille = [dpdx, mu](const sctl::Vector<Real> X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        U.SetZero();
        for (sctl::Long i = 0; i < N; i++) {
            sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
            U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4./mu; // TODO: check direction
        }
        return U;
    };

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
    // std::cout << "size inner = " << size_inner << ", size outer = " << size_outer << std::endl;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }

    sctl::Vector<Real> X1temp, X2temp;
    sinusoid.GetInnerCoord(&X1temp, nullptr);
    sinusoid.GetOuterCoord(&X2temp, nullptr);
    std::cout << "size of X1 temp: " << X1temp.Dim() << ", size of X2 temp: " << X2temp.Dim() << std::endl;
    sinusoid.WriteVTK("../vis/sine_channel_inner", "../vis/sine_channel_outer", X1temp, X2temp, comm);

    // Create lambda function
    const auto BIO = [&LPO, D_scal, &NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real> sigma) {
        U->SetZero();
        LPO.ComputePotential(*U, sigma);
        if (D_scal && U->Dim() == sigma.Dim()) (*U) += 0.5*sigma*D_scal*NormalOrient;
    };

    // Setup gmres
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    // Utot = Ubg + Uwall = 0 on wall. 
    sctl::Vector<Real> vslip = sinusoid.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in_x1, &R_out_x1, dpdx, mu, &sigma, &BIO, &bg_poiseuille](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out_x1 - R_in_x1) + R_in_x1;
        sctl::Vector<Real> r_wts = wts * (R_out_x1 - R_in_x1);

        sctl::Vector<Real> theta_nodes(theta_ord);
        const Real theta_wt = 2*sctl::const_pi<Real>()/theta_ord;
        for (int i=0; i<theta_ord; i++) {
            theta_nodes[i] = theta_wt * i;
        }

        // form mesh grid of x-y points, same-r-all-thetas order.
        sctl::Vector<Real> Xtrg(3 * r_ord * theta_ord);
        for (int i=0; i<r_ord; i++) {
            for (int j=0; j<theta_ord; j++) {
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
        Utrg += bg_poiseuille(Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        // Real Qsa = 0.;
        for (int i=0; i<r_ord; i++) {
            for (int j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * r_nodes[i] * theta_wt * r_wts[i];
                // Qsa += r_nodes[i] * theta_wt * r_wts[i];
            }
        }
        // std::cout << "Debug Qintegral: surface area calculated to be " << Qsa << ", exact is " << sctl::const_pi<Real>() * (0.336*0.336 - 0.13*0.13) << std::endl; // hardcoded to check downstream SA specifically for channel2 geom.

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
    for (int i=0; i<r1_trg.Dim(); i++) {
        r1_trg[i] = getRin(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> r2_trg(ElemOrder_trg*Nelem_trg);
    for (int i=0; i<r2_trg.Dim(); i++) {
        r2_trg[i] = getRout(Xc_trg[i*3+0]);
    }
    sctl::Vector<Real> drdx_trg(ElemOrder_trg*Nelem_trg); // no outward velocity on wall.
    drdx_trg = 0.;
    sctl::Vector<Real> r3 = r1_trg + 1./3. * (r2_trg - r1_trg);
    sctl::Vector<Real> r4 = r1_trg + 2./3. * (r2_trg - r1_trg);
    Annular<Real> channel_trg(Xc_trg,Xc_trg,r3,r4);
    channel_trg.Setup(Nelem_trg, ElemOrder_trg, FourierOrder_trg,drdx_trg,drdx_trg);
    sctl::Vector<Real> X_inner, X_outer;
    channel_trg.GetInnerCoord(&X_inner, nullptr);
    channel_trg.GetOuterCoord(&X_outer, nullptr);
    sctl::Vector<Real> U_inner, U_outer;
    LPO.SetTargetCoord(X_inner);
    BIO(&U_inner, sigma);
    LPO.SetTargetCoord(X_outer);
    BIO(&U_outer, sigma);
    U_inner += bg_poiseuille(X_inner);
    U_outer += bg_poiseuille(X_outer);

    // write to file or read and compare for error. 
    // std::string filename_in = "../out/SelfConv/Concentric_sin_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_U_exact_"+std::to_string(comm.Rank())+"_inner.txt";
    // std::string filename_out = "../out/SelfConv/Concentric_sin_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_U_exact_"+std::to_string(comm.Rank())+"_outer.txt";
    // std::string filename_q = "../out/SelfConv/Concentric_sin_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_Q_exact_"+std::to_string(comm.Rank())+".txt";
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
            // std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
            // std::cout<<"Max relative error = "<< std::setprecision(10) << err_all[0] / u_all[0] << std::endl;
            // std::cout << "Relative error of flux Q = " << std::setprecision(10) << Qerr << std::endl;
        }
    }

}

// */


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

        // concentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm); 

        // eccentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm);

        if (Nelem==4 && FourierOrder == 32) {
            concentric_sine_selfconv<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 1);
        } else {
            concentric_sine_selfconv<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, 0);
        }
        
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

