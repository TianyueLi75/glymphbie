#include <GlymphBIE.hpp>

#include "unit_test_framework.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

// Test 1: Use known Chebyshev quadrature points to test GetCenterLine() function.
template <class Real>
TEST_CASE(test_centerline)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = static_cast<Real>(0.0380602);
    Xc[3] = static_cast<Real>(0.3086583);
    Xc[6] = static_cast<Real>(0.6913417);
    Xc[9] = static_cast<Real>(0.9619398);
    // Cheb from code
    sctl::Vector<Real> Xc_code = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    SCTL_ASSERT(Xc.Dim() == Xc_code.Dim());
    for (sctl::Long ind=0; ind<Xc.Dim(); ind++) {
        ASSERT_NEAR(Xc[ind], Xc_code[ind], 1e-6);
    }
}

template <class Real>
sctl::Vector<Real> read_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    sctl::Vector<Real> data;
    Real x, y, z;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        if (ss >> x >> y >> z) {
            data.PushBack(x);
            data.PushBack(y);
            data.PushBack(z);
        }
    }
    return data;
}

// Test 2: Straight, concentric channels of radii 0.3 and 0.5, along line (x,0.5,0.5), 1 panel and 4 Fourier modes, order = 4. Assume x node values are already at Cheb quadrature points.
// Test setup and accessor
template <class Real>
TEST_CASE(straight_getnodes)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = static_cast<Real>(0.0380602);
    Xc[3] = static_cast<Real>(0.3086583);
    Xc[6] = static_cast<Real>(0.6913417);
    Xc[9] = static_cast<Real>(0.9619398);
    sctl::Vector<Real> r1(ElemOrder*Nelem);
    r1 = 0.3;
    sctl::Vector<Real> r2(ElemOrder*Nelem);
    r2 = 0.5;
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt); 
    
    Real rad_out = straight.GetMinRadius();
    ASSERT_NEAR(rad_out, r1[0], 1e-6);
    
    sctl::Vector<Real> X_inner, X_outer, X_all, Xn_inner, Xn_outer, Xn_all;
    straight.GetInnerCoord(&X_inner, &Xn_inner);
    straight.GetOuterCoord(&X_outer, &Xn_outer);
    straight.GetNodeCoord(&X_all, &Xn_all);
    sctl::Vector<Real> X_file_inner = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_r=0pt3_Np=1_Cheb=4_Nf=4.txt"); // file combines nodes then normal
    sctl::Vector<Real> X_file_outer = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_r=0pt5_Np=1_Cheb=4_Nf=4.txt");
    sctl::Vector<Real> X_file_all = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_annular_Np=1_Cheb=4_Nf=4.txt"); 

    SCTL_ASSERT(X_file_inner.Dim() == X_inner.Dim() + Xn_inner.Dim());
    SCTL_ASSERT(X_file_outer.Dim() == X_outer.Dim() + Xn_outer.Dim());
    SCTL_ASSERT(X_file_all.Dim() == X_all.Dim() + Xn_all.Dim());

    SCTL_ASSERT(X_all.Dim() == X_inner.Dim() + X_outer.Dim());
    SCTL_ASSERT(Xn_all.Dim() == Xn_inner.Dim() + Xn_outer.Dim());
    
    // std::cout << "file inner" << std::endl;
    for (int i=0; i<X_file_inner.Dim(); i++) {
        if (i < X_inner.Dim()) {
            ASSERT_NEAR(X_inner[i], X_file_inner[i], 1e-6);
        } else {
            sctl::Long in = i - X_inner.Dim();
            // std::cout << "Xn value is " << Xn_inner[in] << ", from file is " << X_file_inner[i] << std::endl;
            ASSERT_NEAR(Xn_inner[in], X_file_inner[i], 1e-6);
        }
    }
    // std::cout << "file outer" << std::endl;
    for (int i=0; i<X_file_outer.Dim(); i++) {
        if (i < X_outer.Dim()) {
            // std::cout << "X value is " << X_outer[i] << ", from file is " << X_file_outer[i] << std::endl;
            ASSERT_NEAR(X_outer[i], X_file_outer[i], 1e-6);
        } else {
            sctl::Long in = i - X_outer.Dim();
            // std::cout << "Xn value is " << Xn_outer[in] << ", from file is " << X_file_outer[i] << std::endl;
            ASSERT_NEAR(Xn_outer[in], X_file_outer[i], 1e-6);
        }
    }
    // std::cout << "file all" << std::endl;
    for (int i=0; i<X_file_all.Dim(); i++) {
        if (i < X_all.Dim()) {
            // std::cout << "X value is " << X_all[i] << ", from file is " << X_file_all[i] << std::endl;
            ASSERT_NEAR(X_all[i], X_file_all[i], 1e-6);
        } else {
            sctl::Long in = i - X_all.Dim();
            // std::cout << "X value is " << Xn_all[in] << ", from file is " << X_file_all[i] << std::endl;
            ASSERT_NEAR(Xn_all[in], X_file_all[i], 1e-6);
        }
    }
}

// Test 3.1: check interpolation: 1) fixed r, wrong number of nodes for x, check error output and new x. 2) varying r, test interpolation with precomputed results
template <class Real>
TEST_CASE(interp1)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;

    // Centerline and radii information that don't match the Cheb nodes
    sctl::Vector<Real> Xc(3);
    Xc = 0.5;
    sctl::Vector<Real> r1(1);
    r1 = 0.3;
    sctl::Vector<Real> r2(1);
    r2 = 0.5;
    sctl::Vector<Real> drdt(1); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    auto [Xc1new, Xc2new, r1new, r2new, drdt1new, drdt2new] = straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);

    Real rad_out = straight.GetMinRadius();
    ASSERT_NEAR(rad_out, r1new[0], 1e-6);

    SCTL_ASSERT(Xc1new.Dim()==3*ElemOrder);
    SCTL_ASSERT(Xc2new.Dim()==3*ElemOrder);
    SCTL_ASSERT(r1new.Dim()==ElemOrder);
    SCTL_ASSERT(r2new.Dim()==ElemOrder);
    SCTL_ASSERT(drdt1new.Dim()==ElemOrder);
    SCTL_ASSERT(drdt2new.Dim()==ElemOrder);
    
    sctl::Vector<Real> X_inner, X_outer, X_all, Xn_inner, Xn_outer, Xn_all;
    straight.GetInnerCoord(&X_inner, &Xn_inner);
    straight.GetOuterCoord(&X_outer, &Xn_outer);
    straight.GetNodeCoord(&X_all, &Xn_all);
    sctl::Vector<Real> X_file_inner = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_r=0pt3_Np=1_Cheb=4_Nf=4.txt"); // file combines nodes then normal
    sctl::Vector<Real> X_file_outer = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_r=0pt5_Np=1_Cheb=4_Nf=4.txt");
    sctl::Vector<Real> X_file_all = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_annular_Np=1_Cheb=4_Nf=4.txt"); 

    SCTL_ASSERT(X_file_inner.Dim() == X_inner.Dim() + Xn_inner.Dim());
    SCTL_ASSERT(X_file_outer.Dim() == X_outer.Dim() + Xn_outer.Dim());
    SCTL_ASSERT(X_file_all.Dim() == X_all.Dim() + Xn_all.Dim());

    SCTL_ASSERT(X_all.Dim() == X_inner.Dim() + X_outer.Dim());
    SCTL_ASSERT(Xn_all.Dim() == Xn_inner.Dim() + Xn_outer.Dim());
    
    // std::cout << "file inner" << std::endl;
    for (int i=0; i<X_file_inner.Dim(); i++) {
        if (i < X_inner.Dim()) {
            ASSERT_NEAR(X_inner[i], X_file_inner[i], 1e-6);
        } else {
            sctl::Long in = i - X_inner.Dim();
            // std::cout << "Xn value is " << Xn_inner[in] << ", from file is " << X_file_inner[i] << std::endl;
            ASSERT_NEAR(Xn_inner[in], X_file_inner[i], 1e-6);
        }
    }
    // std::cout << "file outer" << std::endl;
    for (int i=0; i<X_file_outer.Dim(); i++) {
        if (i < X_outer.Dim()) {
            // std::cout << "X value is " << X_outer[i] << ", from file is " << X_file_outer[i] << std::endl;
            ASSERT_NEAR(X_outer[i], X_file_outer[i], 1e-6);
        } else {
            sctl::Long in = i - X_outer.Dim();
            // std::cout << "Xn value is " << Xn_outer[in] << ", from file is " << X_file_outer[i] << std::endl;
            ASSERT_NEAR(Xn_outer[in], X_file_outer[i], 1e-6);
        }
    }
    // std::cout << "file all" << std::endl;
    for (int i=0; i<X_file_all.Dim(); i++) {
        if (i < X_all.Dim()) {
            // std::cout << "X value is " << X_all[i] << ", from file is " << X_file_all[i] << std::endl;
            ASSERT_NEAR(X_all[i], X_file_all[i], 1e-6);
        } else {
            sctl::Long in = i - X_all.Dim();
            // std::cout << "X value is " << Xn_all[in] << ", from file is " << X_file_all[i] << std::endl;
            ASSERT_NEAR(Xn_all[in], X_file_all[i], 1e-6);
        }
    }

}

// Test 3.2: check interpolation: 1) fixed r, wrong number of nodes for x, check error output and new x. 2) varying r, test interpolation with precomputed results
template <class Real>
TEST_CASE(interp2)
{
    // Specify parameters
    sctl::Long ElemOrder = 20;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;

    const auto r1f = [](const sctl::Vector<Real> Xvec) {
        sctl::Vector<Real> rvec(Xvec.Dim()/3);
        for (sctl::Long ind=0; ind<rvec.Dim(); ind++) {
            Real xhere = Xvec[ind*3+0];
            rvec[ind] = 0.5*sctl::sin<Real>(2*sctl::const_pi<Real>()*xhere)+2;
        }
        return rvec;
    };
    const auto r2f = [](const sctl::Vector<Real> Xvec) {
        sctl::Vector<Real> rvec(Xvec.Dim()/3);
        for (sctl::Long ind=0; ind<rvec.Dim(); ind++) {
            Real xhere = Xvec[ind*3+0];
            rvec[ind] = 1.0*sctl::sin<Real>(2*sctl::const_pi<Real>()*Xvec[ind*3+0])+2;
        }
        return rvec;
    };

    // Centerline and radii information that don't match the Cheb nodes
    sctl::Long OldOrder = 16;
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, OldOrder);
    sctl::Vector<Real> r1 = r1f(Xc);
    sctl::Vector<Real> r2 = r2f(Xc);
    sctl::Vector<Real> drdt(1 * OldOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    auto [Xc1new, Xc2new, r1new, r2new, drdt1new, drdt2new] = straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt,true);
    
    sctl::Vector<Real> Xc_code = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r3 = r1f(Xc_code);
    sctl::Vector<Real> r4 = r2f(Xc_code);

    SCTL_ASSERT(Xc1new.Dim()==Xc_code.Dim());
    SCTL_ASSERT(Xc2new.Dim()==Xc_code.Dim());
    SCTL_ASSERT(r1new.Dim()==r3.Dim());
    SCTL_ASSERT(r2new.Dim()==r4.Dim());

    for (int i=0; i<Xc1new.Dim(); i++) {
        ASSERT_NEAR(Xc1new[i], Xc_code[i], 1e-6);
        ASSERT_NEAR(Xc2new[i], Xc_code[i], 1e-6);
    }
    for (int i=0; i<r1new.Dim(); i++) {
        ASSERT_NEAR(r1new[i], r3[i], 1e-6);
    }
    for (int i=0; i<r2new.Dim(); i++) {
        ASSERT_NEAR(r2new[i], r4[i], 1e-6);
    }

}

// Test 3.3: Multiple panel interp.
template <class Real>
TEST_CASE(interp3)
{
    // Specify parameters
    sctl::Long ElemOrder = 16;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;

    const auto r1f = [](const sctl::Vector<Real> Xvec) {
        sctl::Vector<Real> rvec(Xvec.Dim()/3);
        for (sctl::Long ind=0; ind<rvec.Dim(); ind++) {
            Real xhere = Xvec[ind*3+0];
            // rvec[ind] = 0.5*sctl::sin<Real>(2*sctl::const_pi<Real>()*xhere)+2;
            rvec[ind] = 0.5+0.5*xhere;
        }
        return rvec;
    };
    const auto r2f = [](const sctl::Vector<Real> Xvec) {
        sctl::Vector<Real> rvec(Xvec.Dim()/3);
        for (sctl::Long ind=0; ind<rvec.Dim(); ind++) {
            Real xhere = Xvec[ind*3+0];
            // rvec[ind] = 1.0*sctl::sin<Real>(2*sctl::const_pi<Real>()*Xvec[ind*3+0])+2;
            rvec[ind] = 1.0 + 0.5*xhere;
        }
        return rvec;
    };

    // Centerline and radii information that don't match the Cheb nodes
    sctl::Long OldOrder = 12;
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, OldOrder);
    sctl::Vector<Real> r1 = r1f(Xc);
    sctl::Vector<Real> r2 = r2f(Xc);
    sctl::Vector<Real> drdt(1 * OldOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    auto [Xc1new, Xc2new, r1new, r2new, drdt1new, drdt2new] = straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt,true);
    
    sctl::Vector<Real> Xc_code = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r3 = r1f(Xc_code);
    sctl::Vector<Real> r4 = r2f(Xc_code);

    SCTL_ASSERT(Xc1new.Dim()==Xc_code.Dim());
    SCTL_ASSERT(Xc2new.Dim()==Xc_code.Dim());
    SCTL_ASSERT(r1new.Dim()==r3.Dim());
    SCTL_ASSERT(r2new.Dim()==r4.Dim());

    for (int i=0; i<Xc1new.Dim(); i++) {
        ASSERT_NEAR(Xc1new[i], Xc_code[i], 1e-6);
        ASSERT_NEAR(Xc2new[i], Xc_code[i], 1e-6);
    }
    for (int i=0; i<r1new.Dim(); i++) {
        // std::cout << "debug first view: " << std::setprecision(8) << r1new[i] << ", " << r3[i] << std::endl;
        ASSERT_NEAR(r1new[i], r3[i], 1e-6);
    }
    for (int i=0; i<r2new.Dim(); i++) {
        // std::cout << "debug first view, outer: " << std::setprecision(8) << r2new[i] << ", " << r4[i] << std::endl;
        ASSERT_NEAR(r2new[i], r4[i], 1e-6);
    }

}

// Test 4: in_domain check
template <class Real>
TEST_CASE(in_domain)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = static_cast<Real>(0.0380602);
    Xc[3] = static_cast<Real>(0.3086583);
    Xc[6] = static_cast<Real>(0.6913417);
    Xc[9] = static_cast<Real>(0.9619398);
    sctl::Vector<Real> r1(ElemOrder);
    r1 = 0.3;
    sctl::Vector<Real> r2(ElemOrder);
    r2 = 0.5;
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);

    sctl::Vector<Real> Xtrg(3);
    Xtrg = {0.75,0.75,0.75};
    sctl::Vector<sctl::Long> in_dom = straight.InDomain(Xtrg);
    SCTL_ASSERT(in_dom.Dim()==1);
    SCTL_ASSERT(in_dom[0]);
}

// Test 5: Get vslip 
template <class Real>
TEST_CASE(get_vslip) 
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = static_cast<Real>(0.0380602);
    Xc[3] = static_cast<Real>(0.3086583);
    Xc[6] = static_cast<Real>(0.6913417);
    Xc[9] = static_cast<Real>(0.9619398);
    sctl::Vector<Real> r1(ElemOrder);
    r1 = 0.3;
    sctl::Vector<Real> r2(ElemOrder);
    r2 = 0.5;
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    
    // First, this should return all zeros.
    drdt = 0.;
    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);
    sctl::Vector<Real> v1 = straight.GetVslip();
    for (int i=0; i<v1.Dim(); i++) {
        ASSERT_NEAR(v1[i], 0., 1e-6);
    }
    std::cout << "passed first all 0 test." << std::endl;

    //Next, nonzero radial change. FourierOrder = 4 ==> just four directions
    drdt = 0.5;
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt, true);
    sctl::Vector<Real> v2 = straight.GetVslip();
    sctl::Vector<Real> vfile = read_from_file<Real>("../../tests/Annular_tests_files/straight_channel_drdt.txt"); 

    // std::cout << "====================== get vslip values: " << std::endl;
    // for (int i=0; i<v2.Dim(); i++) {
    //     std::cout << v2[i] << std::endl;
    // }
    // std::cout << "====================== 'true' values: " << std::endl;
    // for (int i=0; i<vfile.Dim(); i++) {
    //     std::cout << vfile[i] << std::endl;
    // }

    SCTL_ASSERT(v2.Dim() == vfile.Dim());
    for (int i=0; i<vfile.Dim(); i++) {
        ASSERT_NEAR(vfile[i], v2[i], 1e-6);
    }
    
}

// Test 6: SetupInner should not change Xc and r if already set up and not "force_setup"
template <class Real>
TEST_CASE(setup_booleans)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = static_cast<Real>(0.0380602);
    Xc[3] = static_cast<Real>(0.3086583);
    Xc[6] = static_cast<Real>(0.6913417);
    Xc[9] = static_cast<Real>(0.9619398);
    sctl::Vector<Real> r1(ElemOrder);
    r1 = 0.3;
    sctl::Vector<Real> r2(ElemOrder);
    r2 = 0.5;
    sctl::Vector<Real> drdt(ElemOrder*Nelem); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt); 

    sctl::Vector<Real> Xc_inner_2, Xc_outer_2, r_inner_2, r_outer_2, drdt_inner_2, drdt_outer_2;
    std::tie(Xc_inner_2, r_inner_2, drdt_inner_2) = straight.SetupInner(Nelem, ElemOrder, FourierOrder, drdt);
    std::tie(Xc_outer_2, r_outer_2, drdt_outer_2) = straight.SetupOuter(Nelem, ElemOrder, FourierOrder, drdt);

    SCTL_ASSERT(Xc_inner_2.Dim() == Xc.Dim());
    SCTL_ASSERT(r_inner_2.Dim() == r1.Dim());
    SCTL_ASSERT(r_outer_2.Dim() == r2.Dim());
    SCTL_ASSERT(drdt_inner_2.Dim() == drdt.Dim());

    // Tolerance is really low here since the values should not change at all.
    for (int i=0; i<Xc_inner_2.Dim(); i++) {
        ASSERT_NEAR(Xc_inner_2[i], Xc[i], 1e-8);
        ASSERT_NEAR(Xc_outer_2[i], Xc[i], 1e-8);
    }
    for (int i=0; i<r_inner_2.Dim(); i++) {
        ASSERT_NEAR(r_inner_2[i], r1[i], 1e-8);
        ASSERT_NEAR(r_outer_2[i], r2[i], 1e-8);
    }
    for (int i=0; i<drdt_inner_2.Dim(); i++) {
        ASSERT_NEAR(drdt_inner_2[i], drdt[i], 1e-8);
        ASSERT_NEAR(drdt_outer_2[i], drdt[i], 1e-8);
    }
}

// Test 7: Setter and Getter for centerline nodes and radii; setup booleans test
template <class Real>
TEST_CASE(update_centerline)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    sctl::Long FourierOrder = 4;
    // Make Xc and r
    sctl::Vector<Real> Xc = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(Xc.Dim()/3);
    r1 = 0.3;
    sctl::Vector<Real> r2(Xc.Dim()/3);
    r2 = 0.5;
    sctl::Vector<Real> drdt(Xc.Dim()/3); // no outward velocity on wall.
    drdt = 0.;

    Annular<Real> straight(Xc,Xc,r1,r2);
     sctl::Vector<Real> Xc_inner, Xc_outer, r_inner, r_outer, drdt_inner, drdt_outer;
    std::tie(Xc_inner, Xc_outer, r_inner, r_outer, drdt_inner, drdt_outer) = straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt);
    
    // Update order of centerline by making a new Xc and setting inner centerline to new values.
    sctl::Vector<Real> Xc_upsample = Annular<Real>::GetCenterLine(Nelem*2, ElemOrder);
    sctl::Vector<Real> r3(Xc_upsample.Dim()/3);
    r3 = r1[0];
    sctl::Vector<Real> drdt_upsample(Xc_upsample.Dim()/3);
    drdt_upsample = drdt[0];

    straight.SetInnerXc(Xc_upsample);
    straight.SetInnerR(r3);
    sctl::Vector<Real> Xc_inner_2, Xc_outer_2, r_inner_2, r_outer_2, drdt_inner_2, drdt_outer_2;
    std::tie(Xc_inner_2, r_inner_2, drdt_inner_2) = straight.SetupInner(Nelem*2, ElemOrder, FourierOrder, drdt_upsample); 
    std::tie(Xc_outer_2, r_outer_2, drdt_outer_2) = straight.SetupOuter(Nelem, ElemOrder, FourierOrder, drdt);

    SCTL_ASSERT(Xc_inner_2.Dim() == Xc_upsample.Dim());
    SCTL_ASSERT(Xc_outer_2.Dim() == Xc_outer.Dim());
    SCTL_ASSERT(r_inner_2.Dim() == r3.Dim());
    SCTL_ASSERT(r_outer_2.Dim() == r_outer.Dim());
    SCTL_ASSERT(drdt_inner_2.Dim() == drdt_upsample.Dim());
    SCTL_ASSERT(drdt_outer_2.Dim() == drdt.Dim());

    // Tolerance is really low here since the values should not change at all.
    for (int i=0; i<Xc_inner_2.Dim(); i++) {
        ASSERT_NEAR(Xc_inner_2[i], Xc_upsample[i], 1e-8);
    }
    for (int i=0; i<Xc_outer_2.Dim(); i++) {
        ASSERT_NEAR(Xc_outer_2[i], Xc_outer[i], 1e-8);
    }
    for (int i=0; i<r_inner_2.Dim(); i++) {
        ASSERT_NEAR(r_inner_2[i], r3[i], 1e-8);
    }
    for (int i=0; i<r_outer_2.Dim(); i++) {
        ASSERT_NEAR(r_outer_2[i], r_outer[i], 1e-8);
    }
    for (int i=0; i<drdt_inner_2.Dim(); i++) {
        ASSERT_NEAR(drdt_inner_2[i], drdt_upsample[i], 1e-8);
    }
    for (int i=0; i<drdt_outer_2.Dim(); i++) {
        ASSERT_NEAR(drdt_outer_2[i], drdt_outer[i], 1e-8);
    }

}




// Define test suite
template <class Real>
TEST_SUITE(TestSuite1)
{
    TEST(test_centerline<Real>);
    TEST(straight_getnodes<Real>);
    TEST(interp1<Real>);
    TEST(in_domain<Real>);
    TEST(get_vslip<Real>);
    TEST(interp2<Real>);
    TEST(interp3<Real>);
    TEST(setup_booleans<Real>);
    TEST(update_centerline<Real>);
}

// MPI Tests comparing to serial code.
template <class Real>
void test_centerline(sctl::Comm comm)
{
    // std::cout << "\n Starting test_centerline test." << std::endl;
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

    // std::cout << "Passed test_centerline test." << std::endl;
}

template <class Real>
void straight_getnodes(sctl::Comm comm) 
{
    // std::cout << "\n Starting get_nodes test." << std::endl;
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
    sctl::Vector<Real> drdt(Xc_serial.Dim()/3); // no outward velocity on wall.
    drdt = 0.;
    Annular<Real> straight(Xc_serial,Xc_serial,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt); 

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    // std::cout << "on rank " << comm.Rank() << ", cnt is " << cnt << ", dsp is " << dsp << ", size of Xc_mpi here is " << Xc_mpi.Dim() << std::endl;
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdt_mpi(Xc_mpi.Dim()/3); // no outward velocity on wall.
    drdt_mpi = 0.;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder,drdt_mpi,drdt_mpi);
    
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

    // std::cout << "Passed get_nodes test." << std::endl;
}

template <class Real>
void in_domain(sctl::Comm comm) 
{
    // std::cout << "\n Starting in_domain test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdt_mpi(Xc_mpi.Dim()/3); // no outward velocity on wall.
    drdt_mpi = 0.;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder, drdt_mpi, drdt_mpi);

    sctl::Vector<Real> Xtrg(3);
    Xtrg = {0.75,0.75,0.75};
    sctl::Vector<sctl::Long> in_dom = straight_mpi.InDomain_mpi(Xtrg);
    std::cout << "On rank " << comm.Rank() << "; in domain is " << in_dom[0] << std::endl;

    // Collect in_domain boolean from all process
    sctl::Vector<sctl::Long> in_dom_loc = in_dom;
    sctl::Vector<sctl::Long> in_dom_all(in_dom.Dim());
    in_dom_all = 0;
    comm.Allreduce((sctl::Iterator<Real>) in_dom_loc.begin(), (sctl::Iterator<Real>) in_dom_all.begin(), in_dom.Dim(), sctl::CommOp::MAX); // will be 1 if ==1 on any process; =0 if all ==0.
    in_dom = in_dom_all;

    SCTL_ASSERT(in_dom.Dim()==1);
    SCTL_ASSERT(in_dom[0]);

    // std::cout << "Passed in_domain test." << std::endl;
}

template <class Real>
void get_vslip(sctl::Comm comm) 
{
    // std::cout << "\n Starting get_vslip test." << std::endl;
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 2;
    sctl::Long FourierOrder = 4;

    auto Xc_serial = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    sctl::Vector<Real> r1(Xc_serial.Dim()/3);
    r1 = 0.3;
    sctl::Vector<Real> r2(Xc_serial.Dim()/3);
    r2 = 0.5;
    sctl::Vector<Real> drdt(Xc_serial.Dim()/3); // no outward velocity on wall.
    drdt = 0.5;
    Annular<Real> straight(Xc_serial,Xc_serial,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder,drdt,drdt); 
    sctl::Vector<Real> v2 = straight.GetVslip();

    auto [cnt, dsp, Xc_mpi] = Annular<Real>::GetCenterLine_mpi(Nelem, ElemOrder, comm);
    // std::cout << "on rank " << comm.Rank() << ", cnt is " << cnt << ", dsp is " << dsp << ", size of Xc_mpi here is " << Xc_mpi.Dim() << std::endl;
    sctl::Vector<Real> r1_mpi(Xc_mpi.Dim()/3);
    r1_mpi = 0.3;
    sctl::Vector<Real> r2_mpi(Xc_mpi.Dim()/3);
    r2_mpi = 0.5;
    sctl::Vector<Real> drdt_mpi(Xc_mpi.Dim()/3); 
    drdt_mpi = 0.5;
    Annular<Real> straight_mpi(Xc_mpi,Xc_mpi,r1_mpi,r2_mpi, comm);
    straight_mpi.Setup_mpi(Nelem, ElemOrder, FourierOrder, drdt_mpi, drdt_mpi);
    sctl::Vector<Real> v2_mpi = straight_mpi.GetVslip_mpi();

    sctl::Long offset = dsp*ElemOrder*3*FourierOrder;
    for (sctl::Long ind=0; ind < v2_mpi.Dim(); ind++) {
        SCTL_ASSERT(sctl::fabs(v2_mpi[ind] - v2[offset + ind]) < 1e-6);
    }

    // std::cout << "Passed get_vslip test." << std::endl;
}


int main(int argc, char** argv)
{
    sctl::Comm::MPI_Init(&argc, &argv);
    {
        // Run the unit tests. If a test fails, the program will print failure info
        // and return 1.
        RUN_SUITE(TestSuite1<float>);
        RUN_SUITE(TestSuite1<double>);

        // Test MPI constructions
        sctl::Comm comm = sctl::Comm::World();
        test_centerline<double>(comm);
        straight_getnodes<double>(comm);
        in_domain<double>(comm);
        get_vslip<double>(comm);

        test_centerline<float>(comm);
        straight_getnodes<float>(comm);
        in_domain<float>(comm);
        get_vslip<float>(comm);
    }

    sctl::Comm::MPI_Finalize();
    
    return 0; 
}