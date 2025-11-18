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
    sctl::Vector<Real> r1(ElemOrder);
    r1 = 0.3;
    sctl::Vector<Real> r2(ElemOrder);
    r2 = 0.5;

    Annular<Real> straight(Xc,Xc,r1,r2);
    straight.Setup(Nelem, ElemOrder, FourierOrder); // TODO: there are some std::cout comments here, should also check. E.g. here "cursory check ..." should print. 
    sctl::Vector<Real> X_inner, X_outer, X_all, Xn_inner, Xn_outer, Xn_all;
    straight.GetInnerCoord(&X_inner, &Xn_inner);
    straight.GetOuterCoord(&X_outer, &Xn_outer);
    straight.GetNodeCoord(&X_all, &Xn_all);
    sctl::Vector<Real> X_file_inner = read_from_file<Real>("/home/tianycli/NERS570/F25/glymphbie/tests/Annular_tests_files/straight_channel_r=0pt3_Np=1_Cheb=4_Nf=4.txt"); // file combines nodes then normal
    sctl::Vector<Real> X_file_outer = read_from_file<Real>("/home/tianycli/NERS570/F25/glymphbie/tests/Annular_tests_files/straight_channel_r=0pt5_Np=1_Cheb=4_Nf=4.txt");
    sctl::Vector<Real> X_file_all = read_from_file<Real>("/home/tianycli/NERS570/F25/glymphbie/tests/Annular_tests_files/straight_channel_annular_Np=1_Cheb=4_Nf=4.txt"); 

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
// TODO: more fine grained will test SetupInner_bool, etc. Not now..
// E.g, if calling SetupInner after Setup, will print 'nothing done'

// Test 3: check interpolation: 1) fixed r, wrong number of nodes for x, check error output and new x. 2) varying r, test interpolation with precomputed results


// Test 4: Check SetupInner_bool after SetInnerXc called.





// Define test suite
template <class Real>
TEST_SUITE(TestSuite1)
{
    // TEST(test_centerline<Real>);
    TEST(straight_getnodes<Real>);
}


auto
main() -> int
{
  // Run the unit tests. If a test fails, the program will print failure info
  // and return 1.
  RUN_SUITE(TestSuite1<float>);
//   RUN_SUITE(TestSuite1<double>);
  return 0; 
}