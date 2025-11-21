// Verification for the spatial solver -- Compare with exact solutions, manufactured solutions, and potentially FDM or other published code.
#include <GlymphBIE.hpp>

template <class Real> 
void concentric_poiseuille(const Real R_in, const Real R_out, const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Comm comm)
{
    SCTL_ASSERT(R_out > R_in);

    Real gmres_tol = 1e-10;
    Real tol = 1e-14;

    // Make Annular channel
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(ElemOrder);
    r1 = R_in;
    sctl::Vector<Real> r2(ElemOrder);
    r2 = R_out;
    sctl::Vector<Real> drdx(ElemOrder*Nelem); // no outward velocity on wall.
    drdx = 0.;
    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdx,drdx);
    std::cout << "DEBUG: Xc: "<< std::endl;
    for (int i=0; i<Xc.Dim(); i++) {
        std::cout << Xc[i] << std::endl;
    }

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
                    Real u = -1*coeff * dpdx * (R_out*R_out - r2 - 2 * rmax2 * sctl::log<Real>(rratio));
                    Uvec[starting_ind + 0] = u; // flow always just in the x direction.
                }
            }
        }
        return Uvec;
    };

    // Solve bg poiseuille flow -- dpdx = p(x=1)-p(x=0) = negative for flow right
    const auto bg_poiseuille = [dpdx](const sctl::Vector<Real> X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        U.SetZero();
        for (sctl::Long i = 0; i < N; i++) {
            sctl::Vector<Real> x(3, (sctl::Iterator<Real>) X.begin() + i*3, false);
            U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4; // TODO: check direction
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
    std::cout << "size inner = " << size_inner << ", size outer = " << size_outer << std::endl;
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
    std::cout << "DEBUG, bg poiseuille values, should be equal in magnitude for all nodes on each channel." << std::endl;
    for (int i=0; i<X_annular.Dim()/3; i++) {
        Real magv2 = vbg[i*3+0]*vbg[i*3+0] + vbg[i*3+1]*vbg[i*3+1] + vbg[i*3+2]*vbg[i*3+2];
        std::cout << "vector = " << vbg[i*3+0] << ", " << vbg[i*3+1] << ", " << vbg[i*3+2] << ", magnitude = " << magv2 << std::endl;
        std::cout << "vslip (shoudl be 0): " << vslip[i*3+0] << ", " << vslip[i*3+1] << ", " << vslip[i*3+2] << std::endl;
    }

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compare to exact, report error.
    // Create two cylindrical channels of different radius to test U at different r. (NOTE: Can do sinusoidal channel later)
    sctl::Vector<Real> r3(ElemOrder);
    r3 = R_in + 1./3. * (R_out - R_in);
    sctl::Vector<Real> r4(ElemOrder);
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
    sctl::Vector<Real> Uexact_inner = uexact(X_inner);
    sctl::Vector<Real> Uexact_outer = uexact(X_outer);
    Real max_err_inner = 0.;
    Real max_err_outer = 0.;
    sctl::Vector<Real> Diff_inner = U_inner + bg_poiseuille(X_inner) - Uexact_inner;
    sctl::Vector<Real> Diff_outer = U_outer + bg_poiseuille(X_outer) - Uexact_outer;
    std::cout << " ======================= DEBUG U eval." << std::endl;
    Real mag_exact2 = Uexact_inner[0]*Uexact_inner[0] + Uexact_inner[1]*Uexact_inner[1] + Uexact_inner[2]*Uexact_inner[2];
    std::cout << "uexact = " << Uexact_inner[0] << ", " << Uexact_inner[1] << ", " << Uexact_inner[2] << ", mag = " << mag_exact2 << std::endl;
    Real mag_solve2 = U_inner[0]*U_inner[0] + U_inner[1]*U_inner[1] + U_inner[2]*U_inner[2];
    std::cout << "usolve = " << U_inner[0] << ", " << U_inner[1] << ", " << U_inner[2] << ", mag = " << mag_solve2 << std::endl;
    mag_exact2 = Uexact_outer[0]*Uexact_outer[0] + Uexact_outer[1]*Uexact_outer[1] + Uexact_outer[2]*Uexact_outer[2];
    std::cout << "uexact = " << Uexact_outer[0] << ", " << Uexact_outer[1] << ", " << Uexact_outer[2] << ", mag = " << mag_exact2 << std::endl;
    mag_solve2 = U_outer[0]*U_outer[0] + U_outer[1]*U_outer[1] + U_outer[2]*U_outer[2];
    std::cout << "usolve = " << U_outer[0] << ", " << U_outer[1] << ", " << U_outer[2] << ", mag = " << mag_solve2 << std::endl;
    for (auto e : Diff_inner) max_err_inner = std::max<Real>(max_err_inner, sctl::fabs(e));
    for (auto e : Diff_outer) max_err_outer = std::max<Real>(max_err_outer, sctl::fabs(e));
    std::cout << "Max error at r="<<r3[0]<<": " << max_err_inner << "; max error at r="<<r4[0]<<": " << max_err_outer << std::endl;

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
        double Rin = std::stod(argv[4]); 
        double Rout = std::stod(argv[5]); 
        // Set default physical and problem parameters
        double mu = 1.0;
        double dpdx = -1.;

        concentric_poiseuille<Real>(Rin, Rout, dpdx, mu, Nelem, ElemOrder, FourierOrder, comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

