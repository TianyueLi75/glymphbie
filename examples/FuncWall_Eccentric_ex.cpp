// Verification for the spatial solver -- Compare with exact solutions, manufactured solutions, and potentially FDM or other published code.
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


template <class Real, class InnerFunc, class OuterFunc> 
void flux_sim(const Real dpdx, const Real mu, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, InnerFunc innerFunc,
     OuterFunc outerFunc, Real eccentricity, sctl::Comm comm, const sctl::Long write_ref, const std::string& sim_name)
{
    Real Rin_base = innerFunc.getBaseRadius();
    Real Rout_base = outerFunc.getBaseRadius();
    Real gmres_tol = 1e-11;
    Real tol = 1e-14;
    if (FourierOrder < 16) {
        gmres_tol = 1e-9;
    }

    // const Real R_in = 0.13;
    // const Real R_out = 0.48;

    // // // Channel 1: inner channel wavy, outer fixed
    // // const auto getRin = [&R_in](const Real x) {
    // //     return R_in * (0.3+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    // // };
    // // const auto getRout = [&R_out](const Real x) {
    // //     return R_out;
    // // };

    // // Channel 2: inner channel fixed, outer wavy
    // const auto getRin = [&R_in](const Real x) {
    //     return R_in;
    // };
    // const auto getRout = [&R_out](const Real x) {
    //     return R_out * (0.7+0.2*sctl::sin<Real>(2*sctl::const_pi<Real>()*x));
    // };

    // // Get outlet radius info for flux calculations.
    // Real R_in_x1 = getRin(1.);
    // Real R_out_x1 = getRout(1.);

    // Make Annular channel
    sctl::Vector<Real> Xc_out = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> Xc_in = Xc_out;
    for (sctl::Long i=0; i<Nelem*ElemOrder; i++) {
        Xc_in[i*3+1] += eccentricity; // shift in y-direction
    }
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    // for (int i=0; i<r1.Dim(); i++) {
    //     r1[i] = getRin(Xc[i*3+0]);
    // }
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    // for (int i=0; i<r2.Dim(); i++) {
    //     r2[i] = getRout(Xc[i*3+0]);
    // }
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    // drdt = 0.;
    FuncWall<Real, InnerFunc, OuterFunc> wall(0.01, Xc_in, Xc_out, r1, r2, innerFunc, outerFunc);
    wall.initialize(); // initialize wall geometry before setting up annular geometry
    Annular<Real> annular(wall.centerCoordsIn(),wall.centerCoordsOut(),wall.radiusIn(),wall.radiusOut());
    annular.Setup(Nelem, ElemOrder, FourierOrder,wall.rdotIn(),wall.rdotOut());

    // Make StokesBIO
    // Real min_rad = straight.GetMinRadius();
    // Real S_scal = 1./(2.*min_rad*sctl::log<Real>(1./min_rad)); // TODO: stable computation? check csbq code.
    Real S_scal = 1.0;
    Real D_scal = 1.0;
    StokesBIO<Real> LPO(S_scal, D_scal, comm);
    LPO.AddElemList(annular.GetInnerElemList(), "inner");
    LPO.AddElemList(annular.GetOuterElemList(), "outer");
    LPO.SetAccuracy(tol);
    sctl::Vector<Real> X_annular;
    annular.GetNodeCoord(&X_annular, nullptr);
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

    sctl::Vector<Real> X1temp, X2temp;
    annular.GetInnerCoord(&X1temp, nullptr);
    annular.GetOuterCoord(&X2temp, nullptr);
    std::cout << "size of X1 temp: " << X1temp.Dim() << ", size of X2 temp: " << X2temp.Dim() << std::endl;
    annular.WriteVTK("../vis/"+ sim_name+"_inner", "../vis/" + sim_name+"_outer", X1temp, X2temp, comm);

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
    sctl::Vector<Real> vslip = annular.GetVslip(); // no EXTRA wall velocity
    sctl::Vector<Real> vbg = bg_poiseuille(dpdx,mu,X_annular);

    // Solve gmres
    solver(&sigma, BIO, vslip-vbg, gmres_tol, -1, false, nullptr, &ksp);

    // Compute Q by doing integral over outlet flow
        const auto getQ = [&LPO, &wall,&Rin_base, dpdx, mu, &sigma, &BIO](const sctl::Long r_ord, const sctl::Long theta_ord) {
        SCTL_ASSERT(dpdx<0.);
        sctl::Long last_idx_in = wall.radiusIn().Dim() - 1;
        sctl::Long last_idx_out = wall.radiusOut().Dim() - 1;
        Real R_in = wall.radiusIn()[last_idx_in];
        Real R_out = wall.radiusOut()[last_idx_out];
        Real Xc_y_in = wall.centerCoordsIn()[last_idx_in*3 + 1];
        Real Xc_z_in = wall.centerCoordsIn()[last_idx_in*3 + 2];
        Real Xc_y_out = wall.centerCoordsOut()[last_idx_out*3 + 1];
        Real Xc_z_out = wall.centerCoordsOut()[last_idx_out*3 + 2];

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
        Utrg += bg_poiseuille(dpdx,mu,Xtrg);

        // Quadrature to integrate
        Real Q = 0.;
        for (int i=0; i<theta_ord; i++) {
            for (int j=0; j<r_ord; j++) {
                sctl::Long node_ind = i*r_ord + j;
                Real u = Utrg[3*node_ind + 0];
                Q += u * theta_wt * r_wts_trg[node_ind] * r_nds_trg[node_ind];
            }
        }

        // Real Q2 = Q * mu / (-dpdx) / R_in / R_in / R_in / R_in; // Q in paper = mu/-dp/R_in^4 * int u_code dA should this use R_in_base instead?
        Real Q2 = Q * mu / (-dpdx) / Rin_base / Rin_base / Rin_base / Rin_base;
        return Q2;

    };

    Real Q = getQ(45,64);
    Real eratio = eccentricity / (Rout_base- Rin_base);
    std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder <<  "Eccentricity = " << eratio << "; Flux Q computed = " << Q << " ;" << sim_name << std::endl; // should be the same for all processes
    if (write_ref==1) {
        if (comm.Rank() == 0) { // Only master process should write the final result
            std::string output_dir = "../out/eccentricity_sweep/"; // Ensure this directory exists
            std::string filename = output_dir + sim_name + "_Q_result.txt";
            
            std::ofstream outfile(filename);
            if (outfile.is_open()) {
                outfile << "Dynamics: " << sim_name << "\n";
                outfile << "Eccentricity Ratio: " << eccentricity / (Rout_base- Rin_base) << "\n";
                outfile << "Nelem: " << Nelem << ", ElemOrder: " << ElemOrder << ", FourierOrder: " << FourierOrder << "\n";
                outfile << "Computed Q: " << std::fixed << std::setprecision(10) << Q << "\n";
                outfile.close();
                std::cout << "Wrote Q result to: " << filename << std::endl;
            } else {
                std::cerr << "Error: Could not open output file: " << filename << std::endl;
            }
        }
    }
}
    // // Set target points in between channels 
    // sctl::Long Nelem_trg = 1;
    // sctl::Long ElemOrder_trg = 10;
    // sctl::Long FourierOrder_trg = 16;
    // sctl::Vector<Real> Xc_trg = Annular<Real>::GetCenterLine(Nelem_trg, ElemOrder_trg);
    // sctl::Vector<Real> r1_trg(ElemOrder_trg*Nelem_trg);
    // for (int i=0; i<r1_trg.Dim(); i++) {
    //     r1_trg[i] = getRin(Xc_trg[i*3+0]);
    // }
    // sctl::Vector<Real> r2_trg(ElemOrder_trg*Nelem_trg);
    // for (int i=0; i<r2_trg.Dim(); i++) {
    //     r2_trg[i] = getRout(Xc_trg[i*3+0]);
    // }
    // sctl::Vector<Real> drdt_trg(ElemOrder_trg*Nelem_trg); // no outward velocity on wall.
    // drdt_trg = 0.;
    // sctl::Vector<Real> r3 = r1_trg + 1./3. * (r2_trg - r1_trg);
    // sctl::Vector<Real> r4 = r1_trg + 2./3. * (r2_trg - r1_trg);
    // Annular<Real> channel_trg(Xc_trg,Xc_trg,r3,r4);
    // channel_trg.Setup(Nelem_trg, ElemOrder_trg, FourierOrder_trg,drdt_trg,drdt_trg);
    // sctl::Vector<Real> X_inner, X_outer;
    // channel_trg.GetInnerCoord(&X_inner, nullptr);
    // channel_trg.GetOuterCoord(&X_outer, nullptr);
    // sctl::Vector<Real> U_inner, U_outer;
    // LPO.SetTargetCoord(X_inner);
    // BIO(&U_inner, sigma);
    // LPO.SetTargetCoord(X_outer);
    // BIO(&U_outer, sigma);
    // U_inner += bg_poiseuille(dpdx,mu,X_inner);
    // U_outer += bg_poiseuille(dpdx,mu,X_outer);

    // // write to file or read and compare for error. 
    // std::string filename_in = "../out/SelfConv/Concentric_sin_4_32_U_exact_"+std::to_string(comm.Rank())+"_inner.txt";
    // std::string filename_out = "../out/SelfConv/Concentric_sin_4_32_U_exact_"+std::to_string(comm.Rank())+"_outer.txt";
    // std::string filename_q = "../out/SelfConv/Concentric_sin_4_32_Q_exact_"+std::to_string(comm.Rank())+".txt";
    // if (write_ref==1) {
    //     U_inner.Write(filename_in.c_str());
    //     U_outer.Write(filename_out.c_str());
    //     std::ofstream outfile(filename_q.c_str());
    //     outfile << Q << std::endl;
    //     outfile.close();
    //     // visualization
    //     channel_trg.WriteVTK("../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_inner","../vis/Sine_"+std::to_string(Nelem)+"_"+std::to_string(FourierOrder)+"_outer", U_inner, U_outer, comm);
    // } else {
    //     sctl::Vector<Real> U_ref_in, U_ref_out;
    //     U_ref_in.Read(filename_in.c_str());
    //     U_ref_out.Read(filename_out.c_str());
    //     Real Q_ref;
    //     std::ifstream infile(filename_q.c_str());
    //     infile >> Q_ref;
    //     infile.close();

    //     // Combine error and max_u values on inner and outer channels
    //     sctl::Vector<Real> err = U_inner - U_ref_in;
    //     Real max_err = 0;
    //     for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    //     Real max_u = 0.;
    //     for (const auto e : U_ref_in) max_u = std::max<Real>(max_u, sctl::fabs(e));
    //     err = U_outer - U_ref_out;
    //     for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    //     for (const auto e : U_ref_out) max_u = std::max<Real>(max_u, sctl::fabs(e));

    //     Real Qerr = fabs(Q_ref - Q) / fabs(Q_ref);

    //     sctl::Vector<Real> err_loc(1);
    //     err_loc[0] = max_err;
    //     sctl::Vector<Real> err_all(1);
    //     err_all[0] = 0;
    //     comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);
        
    //     sctl::Vector<Real> u_loc(1);
    //     u_loc[0] = max_u;
    //     sctl::Vector<Real> u_all(1);
    //     u_all[0] = 0.;
    //     comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

    //     if (!comm.Rank()) {
    //         std::cout << "Nelem = " << Nelem << ", Fourier order = " << FourierOrder << "; max error = " << std::setprecision(10) << err_all[0] << "; max relative error = " << err_all[0] / u_all[0] << "; relative error of flux Q is " << Qerr << std::endl;
    //     }
    // }



template <class Real>
struct ConstantFunctor {
    Real _base_radius;

    ConstantFunctor(Real r) : _base_radius(r) {}

    void operator()(sctl::Vector<Real>& radius) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            radius[i] = _base_radius; 
        }
    }

    Real getBaseRadius() const {
        return _base_radius;
    }
};
template <class Real>
struct SpatialSinusoidalFunctor {
    Real _base_radius;
    Real _amplitude;
    Real _frequency;

    SpatialSinusoidalFunctor(Real base_r, Real amp, Real freq) : _base_radius(base_r), _amplitude(amp), _frequency(freq) {}

    void operator()(sctl::Vector<Real>& radius, sctl::Vector<Real> coords) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            Real x = coords[i*3]; // Assuming coords are in (x,y,z) format
            radius[i] = _base_radius * sctl::sqrt( 1 + _amplitude * sctl::sin<Real>(_frequency * sctl::const_pi<Real>() * x) );
        }
    }

    Real getBaseRadius() const {
        return _base_radius;
    }
};

int main(int argc, char** argv)
{
    sctl::Comm::MPI_Init(&argc, &argv);
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <Nelem> <ElemOrder> <FourierOrder> <eccentricity> <inner_func> <outer_func>" << std::endl;
        sctl::Comm::MPI_Finalize();
        return -1;
    }

    using Real=double;

    {
        sctl::Comm comm = sctl::Comm::World();
        // Read parameters from input
        long Nelem = std::stol(argv[1]); // number of elements
        long ElemOrder = std::stol(argv[2]); // Cheb order per element
        long FourierOrder = std::stol(argv[3]);  // number of Fourier nodes
        Real eccentricity_ratio = static_cast<Real>(std::stod(argv[4])); // distance between centers  = eccentricity_ratio * (R_out - R_in)
        long inner_function = std::stol(argv[5]); // options for inner wall shape sinusoidal(1)/circular(0)
        long outer_function = std::stol(argv[6]); // options for outer wall shape sinusoidal(1)/circular(0) 


        // Set default physical and problem parameters
        Real mu = 1.0;
        Real dpdx = -1.;

        // Check input validity
        if (inner_function != 0 && inner_function != 1) {
            if (!comm.Rank()) {
                std::cout << "inner_function option not recognized, please use 0 (circular) or 1 (sinusoidal)." << std::endl;
            }
            sctl::Comm::MPI_Finalize();
            return -1;
        }

        if (outer_function != 0 && outer_function != 1) {
            if (!comm.Rank()) {
                std::cout << "outer_function option not recognized, please use 0 (circular) or 1 (sinusoidal)." << std::endl;
            }
            sctl::Comm::MPI_Finalize();
            return -1;
        }

        Real base_radius_in = 0.3;
        Real amplitude_in = 0.03;
        Real frequency_in = 4.0;
        Real base_radius_out = 0.48;
        Real amplitude_out = 0.05;
        Real frequency_out = 4.0;
        // switch to just based on inner radius?
        Real eccentricity = eccentricity_ratio * (base_radius_out - base_radius_in);
        // --- Define the FOUR possible Functors up front (initialized with a dummy/default state) ---
        SpatialSinusoidalFunctor<Real> inner_sin_functor(base_radius_in, amplitude_in, frequency_in); // r = base_radius * sqrt(1 + amplitude * sin(frequency * pi * x))
        ConstantFunctor<Real> inner_const_functor(base_radius_in); // Using base_radius as the value

        SpatialSinusoidalFunctor<Real> outer_sin_functor(base_radius_out, amplitude_out, frequency_out);
        ConstantFunctor<Real> outer_const_functor(base_radius_out); // Using base_radius as the value

        
        // collision if eccentricity is too large e > (R_out_min - R_in_max)
        const Real R_in_max = base_radius_in * sctl::sqrt(1 + amplitude_in);
        const Real R_out_min = base_radius_out * sctl::sqrt(1 - amplitude_out);
                std::string wall_dynamics_type;
        if (inner_function == 1 && outer_function == 1) {
            // Case: Inner Sinusoidal, Outer Sinusoidal
            wall_dynamics_type = "SS";
        } else if (inner_function == 1 && outer_function == 0) {
            // Case: Inner Sinusoidal, Outer Constant
            wall_dynamics_type = "SC";
        } else if (inner_function == 0 && outer_function == 1) {
            // Case: Inner Constant, Outer Sinusoidal
            wall_dynamics_type = "CS";
        } else { // (inner_function == 0 && outer_function == 0)
            // Case: Inner Constant, Outer Constant
            wall_dynamics_type = "CC";
        }

        std::ostringstream filename_ss;
        filename_ss << wall_dynamics_type 
            << "_E" << std::fixed << std::setprecision(3) << eccentricity_ratio 
            << "_N" << Nelem 
            << "_F" << FourierOrder 
            << "_O" << ElemOrder;
        std::string sim_name = filename_ss.str();


        if (inner_function == 1 && outer_function == 1) {
            // Case: Inner Sinusoidal, Outer Sinusoidal
            SCTL_ASSERT(eccentricity < (R_out_min - R_in_max));
            flux_sim<Real, SpatialSinusoidalFunctor<Real>, SpatialSinusoidalFunctor<Real>>(
                dpdx, mu, Nelem, ElemOrder, FourierOrder, 
                inner_sin_functor, outer_sin_functor, eccentricity,  comm, 1, sim_name
            );
        } else if (inner_function == 1 && outer_function == 0) {
            // Case: Inner Sinusoidal, Outer Constant
            SCTL_ASSERT(eccentricity < (base_radius_out - R_in_max));
            flux_sim<Real, SpatialSinusoidalFunctor<Real>, ConstantFunctor<Real>>(
                dpdx, mu, Nelem, ElemOrder, FourierOrder, 
                inner_sin_functor, outer_const_functor, eccentricity, comm, 1, sim_name
            );
        } else if (inner_function == 0 && outer_function == 1) {
            // Case: Inner Constant, Outer Sinusoidal
            SCTL_ASSERT(eccentricity < (R_out_min - base_radius_in));
            flux_sim<Real, ConstantFunctor<Real>, SpatialSinusoidalFunctor<Real>>(
                dpdx, mu, Nelem, ElemOrder, FourierOrder, 
                inner_const_functor, outer_sin_functor, eccentricity, comm, 1, sim_name
            );
        } else { // (inner_function == 0 && outer_function == 0)
            // Case: Inner Constant, Outer Constant
            SCTL_ASSERT(eccentricity < (base_radius_out - base_radius_in));
            flux_sim<Real, ConstantFunctor<Real>, ConstantFunctor<Real>>(
                dpdx, mu, Nelem, ElemOrder, FourierOrder, 
                inner_const_functor, outer_const_functor, eccentricity, comm, 1, sim_name
            );
        }

    }

    sctl::Comm::MPI_Finalize();
    return 0; 
}

