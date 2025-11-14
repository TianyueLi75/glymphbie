#include <GlymphBIE.hpp>

#include "unit_test_framework.hpp"

// Test 1: Use known GL quadrature points to test GetCenterLine() function.
template <class Real>
TEST_CASE(test_centerline)
{
    // Specify parameters
    sctl::Long ElemOrder = 4;
    sctl::Long Nelem = 1;
    // Make Xc and r
    sctl::Vector<Real> Xc(3*ElemOrder);
    Xc = 0.5; // Base value to make all y,z coordinates 0.5.
    Xc[0] = 0.0694318442;
    Xc[3] = 0.3300094782;
    Xc[6] = 0.6699905218;
    Xc[9] = 0.9305681558;
    // GL from code
    sctl::Vector<Real> Xc_code = Annular<Real>::GetCenterLine(Nelem, ElemOrder);
    SCTL_ASSERT(Xc.Dim() == Xc_code.Dim());
    for (sctl::Long ind=0; ind<Xc.Dim(); ind++) {
        ASSERT_NEAR(Xc[ind], Xc_code[ind], 1e-8);
    }
}

// template <class Real>
// sctl::Vector<Real> read_from_file(std::string filename) {
//     //TODO
// }

// // Test 2: Straight, concentric channels of radii 0.3 and 0.5, along line (x,0.5,0.5), 1 panel and 4 Fourier modes, order = 4. Assume x node values are already at GL quadrature points.
// // Test setup and accessor
// template <class Real>
// TEST_CASE(Test2)
// {
//     // Specify parameters
//     sctl::Long ElemOrder = 4;
//     sctl::Long Nelem = 1;
//     sctl::Long FourierOrder = 4;
//     // Make Xc and r
//     sctl::Vector<Real> Xc(3*ElemOrder);
//     Xc = 0.5; // Base value to make all y,z coordinates 0.5.
//     Xc[0] = 0.0694318442;
//     Xc[3] = 0.3300094782;
//     Xc[6] = 0.6699905218;
//     Xc[9] = 0.9305681558;
//     sctl::Vector<Real> r1(ElemOrder);
//     r1 = 0.3;
//     sctl::Vector<Real> r2(ElemOrder);
//     r2 = 0.5;

//     Annular<Real> straight(Xc,Xc,r1,r2);
//     straight.Setup(Nelem, ElemOrder, FourierOrder); // TODO: there are some std::cout comments here, should also check. E.g. here "cursory check ..." should print. 
//     sctl::Vector<Real> X_inner, X_outer, X_all, Xn_all;
//     straight.GetInnerCoord(X_inner);
//     straight.GetOuterCoord(X_outer);
//     straight.GetNodeCoord(X_all, Xn_all);
//     // TODO: compare to file.
//     sctl::Vector<Real> X_file_inner = read_from_file("straight_channel_r=0pt3_Np=1_GL=4_Nf=4.txt");
//     sctl::Vector<Real> X_file_outer = read_from_file("straight_channel_r=0pt5_Np=1_GL=4_Nf=4.txt");
//     sctl::Vector<Real> X_file_all = read_from_file("straight_channel_annular_Np=1_GL=4_Nf=4.txt"); 
    
// }
// TODO: more fine grained will test SetupInner_bool, etc. Not now..
// E.g, if calling SetupInner after Setup, will print 'nothing done'

// Test 3: check interpolation: 1) fixed r, wrong number of nodes for x, check error output and new x. 2) varying r, test interpolation with precomputed results


// Test 4: Check SetupInner_bool after SetInnerXc called.





// Define test suite
template <class Real>
TEST_SUITE(TestSuite1)
{
    TEST(test_centerline<Real>);
}


auto
main() -> int
{
  // Run the unit tests. If a test fails, the program will print failure info
  // and return 1.
  RUN_SUITE(TestSuite1<float>);
  RUN_SUITE(TestSuite1<double>);
  return 0; 
}