#include "GlymphBIE.hpp"
#include "unit_test_framework.hpp"
#include <cmath>
#include <iostream>
#include <iomanip> // For std::setprecision

// ==========================================
// Functor Definitions
// ==========================================

// 1. Constant Functor (2 Arguments: Radius, Time)
template <class Real>
struct ConstantFunctor {
    Real _target_radius;

    ConstantFunctor(Real r) : _target_radius(r) {}

    void operator()(sctl::Vector<Real>& radius, Real time) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            radius[i] = _target_radius; 
        }
    }
};

// 2. Spatial Functor (3 Arguments: Radius, Time, Coords)
//    Radius varies based on the X-coordinate: R(x) = Base + Slope * x
template <class Real>
struct LinearSpatialFunctor {
    Real _base;
    Real _slope;

    LinearSpatialFunctor(Real base, Real slope) : _base(base), _slope(slope) {}

    // This signature triggers the 'if constexpr' true branch in FuncWall
    void operator()(sctl::Vector<Real>& radius, Real time, const sctl::Vector<Real>& coords) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            // Coords are packed: [x0, y0, z0, x1, y1, z1, ...]
            Real x_coord = coords[i * 3]; 
            radius[i] = _base + _slope * x_coord;
        }
    }
};

// ==========================================
// Test Cases
// ==========================================

template <class Real>
TEST_CASE(constant_walls)
{
    std::cout << "\n[TEST] Running Constant Wall Test..." << std::endl;

    // 1. Setup
    const Real dt = 0.1;
    const sctl::Long N = 5; // Keep grid small for cleaner prints
    const Real target_r_in = 1.0;
    const Real target_r_out = 5.0;

    std::cout << "  -> Initialization: Grid N=" << N << ", dt=" << dt << std::endl;

    sctl::Vector<Real> c_out(N * 3), c_in(N * 3);
    sctl::Vector<Real> r_out(N), r_in(N);
    c_out = 0; c_in = 0;
    
    // Initialize with bad values
    for(sctl::Long i=0; i<N; ++i) { r_out[i] = 99.9; r_in[i] = 99.9; }

    ConstantFunctor<Real> inner_logic(target_r_in);
    ConstantFunctor<Real> outer_logic(target_r_out);

    FuncWall<Real, ConstantFunctor<Real>, ConstantFunctor<Real>> wall(
        dt, c_out, c_in, r_out, r_in, inner_logic, outer_logic
    );

    // 2. Run Update
    std::cout << "  -> Calling update() (Step 1)..." << std::endl;
    wall.update();

    // 3. Verify
    std::cout << "  -> Verifying Time: Expected " << dt << ", Got " << wall.getTime() << std::endl;
    ASSERT_NEAR(wall.getTime(), dt, 1e-10);

    std::cout << "  -> Verifying Geometry (Sample check):" << std::endl;
    bool passed = true;
    for (sctl::Long i = 0; i < N; ++i) {
        if (std::abs(wall.radiusIn()[i] - target_r_in) > 1e-10) passed = false;
        if (std::abs(wall.radiusOut()[i] - target_r_out) > 1e-10) passed = false;
    }
    
    if (passed) std::cout << "     [PASS] All radii matched targets." << std::endl;
    else        std::cout << "     [FAIL] Radius mismatch found." << std::endl;

    for (sctl::Long i = 0; i < N; ++i) {
        ASSERT_NEAR(wall.radiusIn()[i], target_r_in, 1e-10);
        ASSERT_NEAR(wall.radiusOut()[i], target_r_out, 1e-10);
    }
}

template <class Real>
TEST_CASE(spatial_walls)
{
    std::cout << "\n[TEST] Running Spatial (3-Arg) Wall Test..." << std::endl;

    // 1. Setup
    const Real dt = 0.1;
    const sctl::Long N = 5; 
    const Real base_r = 2.0;
    const Real slope = 0.5;

    sctl::Vector<Real> c_in(N * 3), c_out(N * 3); // Only filling 'in' for this test
    sctl::Vector<Real> r_in(N), r_out(N);

    // Setup a grid where X coordinates are 0, 1, 2, 3, 4
    std::cout << "  -> Initialization: Setting up grid x = [0, 1, 2, 3, 4]" << std::endl;
    for(sctl::Long i=0; i<N; ++i) {
        c_in[i*3]     = static_cast<Real>(i); // x
        c_in[i*3 + 1] = 0.0;                  // y
        c_in[i*3 + 2] = 0.0;                  // z
        r_in[i] = 0.0; // Reset radius
    }
    
    // Use Constant for Outer, Spatial for Inner
    ConstantFunctor<Real> outer_logic(10.0); // Keep outer far away
    LinearSpatialFunctor<Real> inner_spatial_logic(base_r, slope);

    FuncWall<Real, LinearSpatialFunctor<Real>, ConstantFunctor<Real>> wall(
        dt, c_out, c_in, r_out, r_in, 
        inner_spatial_logic, outer_logic
    );

    // 2. Run Update
    std::cout << "  -> Calling update()..." << std::endl;
    wall.update();

    // 3. Verify: R should be Base + Slope * x
    //    x=[0,1,2,3,4] -> R=[2.0, 2.5, 3.0, 3.5, 4.0]
    std::cout << "  -> Verifying Spatial Radii:" << std::endl;
    for (sctl::Long i = 0; i < N; ++i) {
        Real x = static_cast<Real>(i);
        Real expected = base_r + slope * x;
        Real actual = wall.radiusIn()[i];

        std::cout << "     Index " << i << " (x=" << x << "): Expected " << expected << ", Got " << actual << std::endl;
        ASSERT_NEAR(actual, expected, 1e-10);
    }
}

// ==========================================
// Test Suite & Main
// ==========================================

TEST_SUITE(func_wall_suite)
{
    TEST(constant_walls<double>);
    TEST(spatial_walls<double>);
}

int main(int argc, char** argv) {
    std::cout << "=========================================" << std::endl;
    std::cout << " STARTING FUNCWALL UNIT TESTS" << std::endl;
    std::cout << "=========================================" << std::endl;
    
    RUN_SUITE(func_wall_suite);
    
    std::cout << "\n=========================================" << std::endl;
    std::cout << " ALL TESTS COMPLETED" << std::endl;
    std::cout << "=========================================" << std::endl;
    return 0;
}