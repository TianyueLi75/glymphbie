// Verification for the spatial solver -- 
// Compare with exact solutions for cylinders.
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

/**
 * @brief For a set of concentric cylinders under constant pressure drop (dp/dx), exact solution is Poiseuille flow.
 *        This test sets up the geometry (hardcoded at this point, see @FuncWall_SpatialSolver_ex for interfacing with Wall)
 *        and evaluates flow at select targets to compute the max relative error from this exact solution.
 *        Result is displayed via print-out.
 *
 * @tparam Real 
 * @param dpdx 
 * @param mu 
 * @param Nelem 
 * @param ElemOrder 
 * @param FourierOrder 
 * @param comm 
 */
template <class Real> 
void concentric_poiseuille(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
{
    const Real R_in = 0.3;
    const Real R_out = 0.48;

    // For Fourier order that is too small, GMRES will plateau early, so a larger tolerance will suffice.
    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    // ==========================
    // Make Annular channel
    // ==========================
    // Set centerline x values that are <ElemOrder>-order Chebyshev nodes
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    // Set radius and wall velocity parameters
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    r1 = R_in;
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    r2 = R_out;
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;
    // Create annular object, concentric with different radii
    Annular<Real> straight(Xc,Xc,r1,r2);
    // Set up annular object using <FourierOrder> nodes in azimuthal direction.
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);
    // Get surface nodes on both channels.
    sctl::Vector<Real> X_annular;
    straight.GetNodeCoord(&X_annular, nullptr);

    // ==========================
    // Exact solution
    // ==========================
    // NOTE: This assumes Xvec taken as an intermediate channel with the same parameters for simplicity. 
    // Otherwise an interpolation is needed.
    const auto uexact = [R_in, R_out, dpdx, mu, &Xc, Nelem, ElemOrder, FourierOrder](const sctl::Vector<Real> Xvec) {
        sctl::Long N = Xvec.Dim()/3;
        SCTL_ASSERT(Xvec.Dim() == Nelem * ElemOrder * FourierOrder * 3);
        sctl::Vector<Real> Uvec(Xvec.Dim());
        Uvec.SetZero();
        // compute constants first
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
                    // Get node coordinates at this index
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

    // ==========================
    // StokesBIO
    // ==========================
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    
    // Use a simple D+S operator, which is full-rank.
    Real S_scal = 1.0; 
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    // Add inner and outer channels to LPO
    LPO.AddElemList(straight.GetInnerElemList(), "inner");
    LPO.AddElemList(straight.GetOuterElemList(), "outer");
    // Set tolerance and periodicity for PVFMM
    LPO.SetAccuracy(tol);
    LPO.SetPeriodicity(sctl::Periodicity::X, 1.0);
    // Set target coordinates and source normal
    LPO.SetTargetCoord(X_annular);
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
    sctl::SlenderElemList<Real> elem_inner = straight.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = straight.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

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

    // ==========================
    // GMRES solve
    // ==========================
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    sctl::Vector<Real> vslip = straight.GetVslip(); // wall velocity from dr/dt
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular); // background flow
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // ==========================
    // Error in flow U
    // ==========================
    // For a easy set of target points, use surface nodes on two cylinders in between.
    sctl::Vector<Real> r3(ElemOrder*Nelem);
    r3 = R_in + 1./3. * (R_out - R_in);
    sctl::Vector<Real> r4(ElemOrder*Nelem);
    r4 = R_in + 2./3. * (R_out - R_in);
    // Set up Annular object to grab surface nodes
    Annular<Real> channel_trg(Xc,Xc,r3,r4);
    channel_trg.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);
    sctl::Vector<Real> X_inner, X_outer;
    channel_trg.GetInnerCoord(&X_inner, nullptr);
    channel_trg.GetOuterCoord(&X_outer, nullptr);
    // Evaluate at target
    sctl::Vector<Real> U_inner, U_outer;
    LPO.SetTargetCoord(X_inner);
    BIO(&U_inner, sigma);
    LPO.SetTargetCoord(X_outer);
    BIO(&U_outer, sigma);
    // Add background flow
    U_inner += bg_poiseuille(dpdx,mu,X_inner);
    U_outer += bg_poiseuille(dpdx,mu,X_outer);
    // Get exact flow
    sctl::Vector<Real> Uexact_inner = uexact(X_inner);
    sctl::Vector<Real> Uexact_outer = uexact(X_outer);
    // Gather max error at either radii, average to one relative error.
    Real max_err_inner = 0.;
    Real max_err_outer = 0.;
    sctl::Vector<Real> Diff_inner = U_inner - Uexact_inner;
    sctl::Vector<Real> Diff_outer = U_outer - Uexact_outer;
    Real max_u_inner = 0.;
    Real max_u_outer = 0.;
    for (auto e : U_inner) max_u_inner = std::max<Real>(max_u_inner, sctl::fabs(e));
    for (auto e : U_outer) max_u_outer = std::max<Real>(max_u_outer, sctl::fabs(e));
    for (auto e : Diff_inner) max_err_inner = std::max<Real>(max_err_inner, sctl::fabs(e));
    for (auto e : Diff_outer) max_err_outer = std::max<Real>(max_err_outer, sctl::fabs(e));
    Real avg_max_err = 0.5*(max_err_inner / max_u_inner + max_err_outer / max_u_outer);
    std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Averaged max relative error at two different r: "<< avg_max_err << std::endl;

    // ==========================
    // Error in flux Q
    // ==========================
    const Real alpha = R_out / R_in;
    const Real beta = alpha;
    const Real K = alpha*beta - 1.0;
    const Real Q_ana = sctl::const_pi<Real>()/8. *((K+1.)*(K+1.) - 1.0 - 2.*K*K/sctl::log<Real>(K+1.));
    
    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in, &R_out, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        // Get Chebyshev node and weights for radial direction
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out - R_in) + R_in;
        sctl::Vector<Real> r_wts = wts * (R_out - R_in);
        // Get equidistant nodes for azimuthal direction
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
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
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

template <class Real> 
void concentric_poiseuille_mpi(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
{
    const Real R_in = 0.3;
    const Real R_out = 0.48;

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    // ==========================
    // Make Annular channel
    // ==========================
    auto [cnt, dsp, Xc] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    sctl::Vector<Real> r1(Xc.Dim()/3);
    r1 = R_in;
    sctl::Vector<Real> r2(Xc.Dim()/3);
    r2 = R_out;
    sctl::Vector<Real> drdt(Xc.Dim()/3); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> straight(Xc,Xc,r1,r2, comm);
    straight.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdt,drdt);

    // ==========================
    // Exact solution
    // ==========================
    // NOTE: This assumes Xvec taken as an intermediate channel with the same parameters for simplicity. 
    // Otherwise an interpolation is needed.
    const auto uexact = [R_in, R_out, dpdx, mu, &Xc, ElemOrder, FourierOrder](const sctl::Vector<Real> Xvec) {
        sctl::Long N = Xvec.Dim()/3;
        SCTL_ASSERT(Xvec.Dim() == Xc.Dim() * FourierOrder);
        sctl::Vector<Real> Uvec(Xvec.Dim());
        Uvec.SetZero();
        // compute constants first 
        Real coeff = 1./4./mu;
        // Get rmax
        Real rmax_num = R_out*R_out - R_in*R_in;
        Real Rratio = R_out / R_in;
        Real rmax_denom = 2*sctl::log<Real>(Rratio);
        Real rmax2 = rmax_num / rmax_denom;
        // Compute u = 1/4/mu*dpdx*(R_out^2-r^2-2rmax^2ln*(R_out/r))
        for (sctl::Long i=0; i<Xc.Dim()/3/ElemOrder; i++) {
            std::cout << "debug in uexact, i = " << i << std::endl;
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

    // ==========================
    // StokesBIO
    // ==========================
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    
    // Use a simple D+S operator, which is full-rank.
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
    sctl::Long size_inner = Xc.Dim() * FourierOrder;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }
    // Periodic single-layer nullspace: far-field quadrature weights and total
    // surface area, used to remove the surface-mean density during the solve.
    sctl::SlenderElemList<Real> elem_inner = straight.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = straight.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

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

    // ==========================
    // GMRES solve
    // ==========================
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> ksp;
    sctl::Vector<Real> sigma;
    sctl::Vector<Real> vslip = straight.GetVslip_mpi(); 
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // ==========================
    // Error in flow U
    // ==========================
    // For a easy set of target points, use surface nodes on two cylinders in between.
    sctl::Vector<Real> r3(Xc.Dim()/3);
    r3 = R_in + 1./3. * (R_out - R_in);
    sctl::Vector<Real> r4(Xc.Dim()/3);
    r4 = R_in + 2./3. * (R_out - R_in);
    Annular<Real> channel_trg(Xc,Xc,r3,r4,comm);
    channel_trg.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdt,drdt);
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
    sctl::Vector<Real> Uexact_inner = uexact(X_inner);
    sctl::Vector<Real> Uexact_outer = uexact(X_outer);
    Real max_err_inner = 0.;
    Real max_err_outer = 0.;
    sctl::Vector<Real> Diff_inner = U_inner - Uexact_inner;
    sctl::Vector<Real> Diff_outer = U_outer - Uexact_outer;
    Real max_u_inner = 0.;
    Real max_u_outer = 0.;
    for (auto e : U_inner) max_u_inner = std::max<Real>(max_u_inner, sctl::fabs(e));
    for (auto e : U_outer) max_u_outer = std::max<Real>(max_u_outer, sctl::fabs(e));
    for (auto e : Diff_inner) max_err_inner = std::max<Real>(max_err_inner, sctl::fabs(e));
    for (auto e : Diff_outer) max_err_outer = std::max<Real>(max_err_outer, sctl::fabs(e));
    Real avg_max_err = 0.5*(max_err_inner / max_u_inner + max_err_outer / max_u_outer);
    // MPI Reduce to collect min of all min radii on each process
    sctl::Vector<Real> err_loc(1);
    err_loc[0] = avg_max_err;
    sctl::Vector<Real> err_all(1);
    err_all[0] = 0.;
    comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);

    if (!comm.Rank()) {
        std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Averaged max relative error at two different r: "<< err_all[0]<< std::endl;
    }

    // Check flux calculation against analytical solutions.
    const Real alpha = R_out / R_in;
    const Real beta = alpha;
    const Real K = alpha*beta - 1.0;
    const Real Q_ana = sctl::const_pi<Real>()/8. *((K+1.)*(K+1.) - 1.0 - 2.*K*K/sctl::log<Real>(K+1.));

    // Compute Q by doing integral over outlet flow
    const auto getQ = [&LPO, &R_in, &R_out, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (R_out - R_in) + R_in;
        sctl::Vector<Real> r_wts = wts * (R_out - R_in);

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
        for (sctl::Long i=0; i<r_ord; i++) {
            for (sctl::Long j=0; j<theta_ord; j++) {
                sctl::Long node_ind = i*theta_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * r_nodes[i] * theta_wt * r_wts[i];
            }
        }

        Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(20,64); // For now just compute flux on all processes

    if (!comm.Rank()) {
        std::cout << "Flux difference: Q analytical = " << std::setprecision(10) << Q_ana << ", Q computed = " << Q << "; rel error is " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;
    }

}

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
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> straight(Xc_inner,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);

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
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }
    // Periodic single-layer nullspace: far-field quadrature weights and total
    // surface area, used to remove the surface-mean density during the solve.
    sctl::SlenderElemList<Real> elem_inner = straight.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = straight.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

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
    sctl::Vector<Real> vslip = straight.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    // r1 and r2 differ at different theta's!
    const auto getQ = [&LPO, &R_in, &R_out, &Xc_y_in, &Xc_y_out, &Xc_z_in, &Xc_z_out, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
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

        for (sctl::Long i=0; i<theta_ord; i++) {
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

            for (sctl::Long j=0; j<r_ord; j++) {
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
        // for (sctl::Long i=0; i<theta_ord; i++) {
        //     for (sctl::Long j=0; j<r_ord; j++) {
        //         sctl::Long node_ind = i*r_ord + j;
        //         Q1 += theta_wt * r_wts_trg[node_ind] * r_nds_trg[node_ind];
        //     }
        // }
        // std::cout << "Q1 from quadrature is " << Q1 << "; SA is " << sctl::const_pi<Real>() * (R_out*R_out - R_in*R_in) << std::endl;

        // Compute U at meshgrid:
        sctl::Vector<Real> Utrg;
        LPO.SetTargetCoord(Xtrg);
        BIO(&Utrg, sigma);
        Utrg += bg_poiseuille(dpdx,mu,Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        for (sctl::Long i=0; i<theta_ord; i++) {
            for (sctl::Long j=0; j<r_ord; j++) {
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

// For eccentric channels, use fixed centerlines and radii to prevent intersection.
template <class Real> 
void eccentric_poiseuille_mpi(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
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
    auto [cnt, dsp, Xc] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    sctl::Vector<Real> Xc_inner = Xc;
    for (sctl::Long i=0; i<Xc_inner.Dim()/3; i++) {
        Xc_inner[i*3+1] = Xc_y_in;
        Xc_inner[i*3+2] = Xc_z_in;
    }
    sctl::Vector<Real> r1(Xc.Dim()/3);
    r1 = R_in;
    sctl::Vector<Real> r2(Xc.Dim()/3);
    r2 = R_out;
    sctl::Vector<Real> drdt(Xc.Dim()/3); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> straight(Xc_inner,Xc,r1,r2, comm);
    straight.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdt,drdt);

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
    sctl::Long size_inner = Xc.Dim() * FourierOrder;
    sctl::Long size_outer = X_annular.Dim() - size_inner;
    for (sctl::Long ind=0; ind<size_outer; ind++) {
        // Normal of outer channel element list points outward by default.
        NormalOrient[size_inner + ind] = -1.;
    }
    // Periodic single-layer nullspace: far-field quadrature weights and total
    // surface area, used to remove the surface-mean density during the solve.
    sctl::SlenderElemList<Real> elem_inner = straight.GetInnerElemList();
    sctl::SlenderElemList<Real> elem_outer = straight.GetOuterElemList();
    sctl::Vector<Real> wts_inner, wts_outer;
    Real surface_area;
    GetSurfWtsArea<Real>(elem_inner, elem_outer, wts_inner, wts_outer, surface_area, comm, tol);

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
    sctl::Vector<Real> vslip = straight.GetVslip_mpi(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
    // r1 and r2 differ at different theta's!
    const auto getQ = [&LPO, &R_in, &R_out, &Xc_y_in, &Xc_y_out, &Xc_z_in, &Xc_z_out, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]

        const auto positiveRoot = [](double a, double b, double c) {
            double disc = b*b - 4*a*c;
            if(disc < 0) {
                std::cout << "discreminant in r_max finding is negative, fails." << std::endl;
                return 0.0; // no real solution
            }
            double r1 = (-b + std::sqrt(disc)) / (2.*a);
            double r2 = (-b - std::sqrt(disc)) / (2.*a);
            double rmax = std::max(r1,r2);
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

        for (sctl::Long i=0; i<theta_ord; i++) {
            Real theta_val = theta_wt * i;

            Real cosT = sctl::cos<Real>(theta_val);
            Real sinT = sctl::sin<Real>(theta_val);
            double b = 2.*(d1*cosT + d2*sinT);
            double c = d2sq - R_out*R_out;
            double rho_max = positiveRoot(1.0, b, c);

            Real r_out = rho_max; 

            sctl::Vector<Real> r_nodes = nds * (r_out - R_in) + R_in;
            sctl::Vector<Real> r_wts = wts * (r_out - R_in);

            r_out_trg[i] = r_out;

            for (sctl::Long j=0; j<r_ord; j++) {
                sctl::Long node_ind = i*r_ord + j;
                Xtrg[3*node_ind + 0] = 1.0; // evaluation at x=1.0 only
                Xtrg[3*node_ind + 1] = Xc_y_in + r_nodes[j] * sctl::cos<Real>(theta_val); 
                Xtrg[3*node_ind + 2] = Xc_z_in + r_nodes[j] * sctl::sin<Real>(theta_val); 

                r_wts_trg[node_ind] = r_wts[j];
                r_nds_trg[node_ind] = r_nodes[j];
            }
        }

        // Compute U at meshgrid:
        sctl::Vector<Real> Utrg;
        LPO.SetTargetCoord(Xtrg);
        BIO(&Utrg, sigma);
        Utrg += bg_poiseuille(dpdx,mu,Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        for (sctl::Long i=0; i<theta_ord; i++) {
            for (sctl::Long j=0; j<r_ord; j++) {
                sctl::Long node_ind = i*r_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * theta_wt * r_wts_trg[node_ind] * r_nds_trg[node_ind];
            }
        }

        Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(45,64);
    if (!comm.Rank()) {
        std::cout <<"Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; Flux relative error = " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;
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

        concentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm); 

        eccentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm);
        
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

