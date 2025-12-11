// Verification for the spatial solver -- Compare with exact solutions, manufactured solutions, and potentially FDM or other published code.
#include <GlymphBIE.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

template <class Real>
struct ConstantFunctor {
    Real _target_radius;

    ConstantFunctor(Real r) : _target_radius(r) {}

    void operator()(sctl::Vector<Real>& radius) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            radius[i] = _target_radius; 
        }
    }
};


template <class Real> 
void concentric_poiseuille(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm, ConstantFunctor<Real> rin_functor, ConstantFunctor<Real> rout_functor)
{
    // const Real R_in = 0.3;
    // const Real R_out = 0.48;

    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    // Make Annular channel
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);

    sctl::Vector<Real> r1(ElemOrder*Nelem); // inner radius
    sctl::Vector<Real> r2(ElemOrder*Nelem); // outer radius

    // rin_functor(r1, 0.0); // initialize wall with this Xc
    // rout_functor(r2, 0.0);
    // set constant radius for now not determined by functor
    r1 = 0.3;
    r2 = 0.48;

    const Real dt = .01;
    FuncWall<Real, ConstantFunctor<Real>, ConstantFunctor<Real>> wall(dt, Xc, Xc, r1, r2, rin_functor, rout_functor);



    // initialize wall with this Xc
    // sctl::Vector<Real> r1(ElemOrder*Nelem);
    // r1 = R_in;
    // sctl::Vector<Real> r2(ElemOrder*Nelem);
    // r2 = R_out;
    // sctl::Vector<Real> drdx(ElemOrder*Nelem); // no outward velocity on wall.
    // drdx = 0.;
    Annular<Real> straight(wall.centerCoordsIn(),wall.centerCoordsOut(),wall.radiusIn(),wall.radiusOut());
    straight.Setup(Nelem, ElemOrder, FourierOrder, wall.rdotIn(), wall.rdotOut());
    // std::cout << "DEBUG: Xc: "<< std::endl;
    // for (int i=0; i<Xc.Dim(); i++) {
    //     std::cout << Xc[i] << std::endl;
    // }
    const Real R_in = r1[0];
    const Real R_out = r2[0];
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
    // this dynamics in functor wall update.
    // sctl::Vector<Real> r3(ElemOrder*Nelem);
    // r3 = R_in + 1./3. * (R_out - R_in); = .36 
    // sctl::Vector<Real> r4(ElemOrder*Nelem);
    // r4 = R_in + 2./3. * (R_out - R_in); = .42
    wall.update();
    Annular<Real> channel_trg(wall.centerCoordsIn(),wall.centerCoordsOut(), wall.radiusIn(),wall.radiusOut());
    channel_trg.Setup(Nelem, ElemOrder, FourierOrder, wall.rdotIn(),wall.rdotOut());
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
    const Real alpha = wall.radiusOut()[0] / wall.radiusIn()[0];
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

    // Compute Q by doing integral over outlet 
    Real Rin = wall.radiusIn()[0];
    Real Rout = wall.radiusOut()[0];
    const auto getQ = [&LPO, &Rin, &Rout, dpdx, mu, &sigma, &BIO, &bg_poiseuille](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        
        sctl::Vector<Real> nds, wts;
        sctl::ChebQuadRule<Real>::ComputeNdsWts(&nds, &wts, r_ord); // nds and wts on [0,1]
        sctl::Vector<Real> r_nodes = nds * (Rout - Rin) + Rin;
        sctl::Vector<Real> r_wts = wts * (Rout - Rin);

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

        Real Q2 = Q * mu / (-dpdx) / Rin / Rin / Rin / Rin; // Q in paper = mu/-dp/R_in^4 * int u_code dA

        return Q2;

    };

    Real Q = getQ(20,64);

    std::cout << "Flux difference: Q analytical = " << std::setprecision(10) << Q_ana << ", Q computed = " << Q << "; rel error is " << fabs(Q-Q_ana) / fabs(Q_ana) << std::endl;

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

        ConstantFunctor<Real> rin_functor(0.3);
        ConstantFunctor<Real> rout_functor(0.48);
        concentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm, rin_functor, rout_functor);
        // eccentric_poiseuille<Real>(dpdx, mu, Nelem, ElemOrder, FourierOrder, comm);

        
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

