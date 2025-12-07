#include "GlymphBIE.hpp"
#include "unit_test_framework.hpp"
#include <cmath>
#include <iostream>
#include <iomanip> // For std::setprecision

template <class Real>
void set_centerline_z_linspace(sctl::Vector<Real>& coords,
                               Real x_fixed,
                               Real y_fixed,
                               Real z_start,
                               Real z_end) {
    const sctl::Long dim = coords.Dim();
    SCTL_ASSERT(dim % 3 == 0);
    const sctl::Long N = dim / 3;

    if (N <= 0) return;
    if (N == 1) {
        coords[0] = x_fixed;
        coords[1] = y_fixed;
        coords[2] = z_start; // or (z_start+z_end)/2
        return;
    }

    const Real dz = (z_end - z_start) / (Real)(N - 1);

    for (sctl::Long i = 0; i < N; ++i) {
        Real z = z_start + dz * (Real)i;
        coords[3*i    ] = x_fixed;
        coords[3*i + 1] = y_fixed;
        coords[3*i + 2] = z;
    }
}



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
//    Radius varies based on the z-coordinate: R(z) = Base + Slope * z
template <class Real>
struct LinearSpatialFunctor {
    Real _base;
    Real _slope;

    LinearSpatialFunctor(Real base, Real slope) : _base(base), _slope(slope) {}

    // This signature triggers the 'if constexpr' true branch in FuncWall
    void operator()(sctl::Vector<Real>& radius, Real time, const sctl::Vector<Real>& coords) const {
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            // Coords are packed: [x0, y0, z0, x1, y1, z1, ...]
            Real z_coord = coords[i * 3+2]; 
            radius[i] = _base + _slope * z_coord;
        }
    }
};

// 3. Spatiotemporal Functor (4 Arguments: Radius, Rdot, Time, Coords)
//    Radius varies based on the z-coordinate: R(z) = Base + Slope * z + Amplitude * sin(omega * t)
template <class Real>
struct SpatiotemporalFunctor {
    Real _base;
    Real _slope;
    Real _amplitude;
    Real _omega;
    SpatiotemporalFunctor(Real base, Real slope, Real amplitude, Real omega) 
        : _base(base), _slope(slope), _amplitude(amplitude), _omega(omega) {}
    // This signature triggers the 'if constexpr' true branch in FuncWall
    void operator()(sctl::Vector<Real>& radius, sctl::Vector<Real>& rdot, Real time, const sctl::Vector<Real>& coords) const {
        Real sin_term = std::sin(_omega * time);
        Real cos_term = std::cos(_omega * time);
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            // Coords are packed: [x0, y0, z0, x1, y1, z1, ...]
            Real z_coord = coords[i * 3+2]; 
            radius[i] = _base + _slope * z_coord + _amplitude * sin_term;
            rdot[i] = _amplitude * _omega * cos_term;
        }
    }
};

// 4. r^2 spatiotemporal Functor (4 Arguments: Radius, Rdot, Time, Coords)
//    Radius varies based on the z-coordinate: r(z,t) = R0*sqrt(1 + Amplitude * cos( 2*pi * n* z/L) * cos(omega * t))
//   Rdot is computed accordingly dr/dt = - (R0 * Amplitude * omega * cos(2*pi*n*z/L) * sin(omega*t)) / (2*sqrt(1 + Amplitude * cos(2*pi*n*z/L) * cos(omega*t)))
template <class Real>
struct R2SpatiotemporalFunctor {
    Real _R0;
    Real _amplitude;
    Real _L;
    Real _omega;
    sctl::Long _n;
    R2SpatiotemporalFunctor(Real R0, Real amplitude, Real L, sctl::Long n, Real omega) 
        : _R0(R0), _amplitude(amplitude), _L(L), _n(n), _omega(omega) {}
    // This signature triggers the 'if constexpr' true branch in FuncWall
    void operator()(sctl::Vector<Real>& radius, sctl::Vector<Real>& rdot, Real time, const sctl::Vector<Real>& coords) const {
        Real cos_omega_t = std::cos(_omega * time);
        Real sin_omega_t = std::sin(_omega * time);
        for (sctl::Long i = 0; i < radius.Dim(); ++i) {
            // Coords are packed: [x0, y0, z0, x1, y1, z1, ...]
            Real z_coord = coords[i * 3+2]; 
            Real cos_term = std::cos(2.0 * sctl::const_pi<Real>() * _n * z_coord / _L);
            Real radicand = 1.0 + _amplitude * cos_term * cos_omega_t;
            radius[i] = _R0 * std::sqrt(radicand);
            // Derivative calculation
            if (radicand > 0.0) {
                rdot[i] = - (_R0 * _amplitude * _omega * cos_term * sin_omega_t) / (2.0 * std::sqrt(radicand));
            } else {
                rdot[i] = 0.0; // Avoid NaN, though this indicates an issue
            }
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
    // Set centerlines (z from 0 to 1)
    set_centerline_z_linspace(c_in, 0.0, 0.0, 0.0, 1.0);
    set_centerline_z_linspace(c_out, 0.0, 0.0, 0.0, 1.0);
    
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

    sctl::Vector<Real> c_in(N * 3), c_out(N * 3); 
    sctl::Vector<Real> r_in(N), r_out(N);

    set_centerline_z_linspace(c_in, 0.0, 0.0, 0.0, 1.0);
    set_centerline_z_linspace(c_out, 0.0, 0.0, 0.0, 1.0);
    r_in = 1.0; r_out = 10.0; // Set radii
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

    // 3. Verify: R should be Base + Slope * z
    std::cout << "  -> Verifying Spatial Radii:" << std::endl;
    for (sctl::Long i = 0; i < N; ++i) {
        Real z = c_in[i * 3 + 2];
        Real expected = base_r + slope * z;
        Real actual = wall.radiusIn()[i];

        std::cout << "     Index " << i << " (z=" << z << "): Expected " << expected << ", Got " << actual << std::endl;
        ASSERT_NEAR(actual, expected, 1e-10);
    }
}

template <class Real>
TEST_CASE(spatiotemporal_walls)
{
    std::cout << "\n[TEST] Running Spatiotemporal (4-Arg) Wall Test..." << std::endl;
    // 1. Setup
    const Real dt = 0.1;
    const sctl::Long N = 5; 
    const Real base_r = 1.0;
    const Real slope = 0.2;
    const Real amplitude = 0.5;
    const Real omega = 2.0 * M_PI; // 1 cycle per unit time 
    sctl::Vector<Real> c_in(N * 3), c_out(N * 3);
    sctl::Vector<Real> r_in(N), r_out(N);   

    // Setup a grid where Z coordinates are linearly spaced [0,1]
    std::cout << "  -> Initialization: Setting up grid z = linspace[0, 1]" << std::endl;
    set_centerline_z_linspace(c_in, 0.0, 0.0, 0.0, 1.0);
    set_centerline_z_linspace(c_out, 0.0, 0.0, 0.0, 1.0);

    r_in = 1.0; r_out = 10.0; // Set radii

    // Use Spatiotemporal for Inner, Constant for Outer
    SpatiotemporalFunctor<Real> inner_func(base_r, slope, amplitude, omega);
    ConstantFunctor<Real> outer_func(10.0); // Keep outer far away
    FuncWall<Real, SpatiotemporalFunctor<Real>, ConstantFunctor<Real>> wall(
        dt, c_out, c_in, r_out, r_in, 
        inner_func, outer_func
    );
    // 2. Run Update
    std::cout << "  -> Calling update()..." << std::endl;
    Real time = wall.getTime();
    wall.update();
    // 3. Verify
    std::cout << "  -> Verifying Spatiotemporal Radii and Velocities:" << std::endl;

    Real sin_term = std::sin(omega * time);
    Real cos_term = std::cos(omega * time);
    for (sctl::Long i = 0; i < N; ++i) {
        Real z = c_in[i * 3 + 2];
        Real expected_radius = base_r + slope * z + amplitude * sin_term;
        Real expected_rdot = amplitude * omega * cos_term;

        Real actual_radius = wall.radiusIn()[i];
        Real actual_rdot = wall.rdotIn()[i];

        std::cout << std::setprecision(10);
        std::cout << "     Index " << i << " (z=" << z << "): "
                  << "Expected R=" << expected_radius << ", Got R=" << actual_radius
                  << "  |  Expected Rdot=" << expected_rdot << ", Got Rdot=" << actual_rdot
                  << std::endl;

        ASSERT_NEAR(actual_radius, expected_radius, 1e-10);
        ASSERT_NEAR(actual_rdot, expected_rdot, 1e-10);
    }
}
template <class Real>
TEST_CASE(r2_spatiotemporal_walls)
{
    std::cout << "\n[TEST] Running r^2 Spatiotemporal (4-Arg) Wall Test..." << std::endl;
    // 1. Setup
    const Real dt = 0.1;
    const sctl::Long N = 5; 
    const Real R0 = 1.0;
    const Real amplitude = 0.3;
    const Real L = 1.0;
    const sctl::Long n = 2;
    const Real omega = 2.0 * M_PI; // 1 cycle per unit time 
    sctl::Vector<Real> c_in(N * 3), c_out(N * 3);
    sctl::Vector<Real> r_in(N), r_out(N);   

    // Setup a grid where Z coordinates are linearly spaced [0,1]
    std::cout << "  -> Initialization: Setting up grid z = linspace[0, 1]" << std::endl;
    set_centerline_z_linspace(c_in, 0.0, 0.0, 0.0, 1.0);
    set_centerline_z_linspace(c_out, 0.0, 0.0, 0.0, 1.0);

    r_in = 1.0; r_out = 10.0; // Set radii

    // Use R2Spatiotemporal for Inner, Constant for Outer
    R2SpatiotemporalFunctor<Real> inner_func(R0, amplitude, L, n, omega);
    ConstantFunctor<Real> outer_func(10.0); // Keep outer far away
    FuncWall<Real, R2SpatiotemporalFunctor<Real>, ConstantFunctor<Real>> wall(
        dt, c_out, c_in, r_out, r_in, 
        inner_func, outer_func
    );
    // 2. Run Update
    std::cout << "  -> Calling update()..." << std::endl;

    for (sctl::Long step = 0; step < 3; ++step) {
        wall.update();
    }
    Real time = wall.getTime() - dt; 
    // 3. Verify
    std::cout << "  -> Verifying r^2 Spatiotemporal Radii and Velocities:" << std::endl;

    Real cos_omega_t = std::cos(omega * time);
    Real sin_omega_t = std::sin(omega * time);
    for (sctl::Long i = 0; i < N; ++i) {
        Real z = c_in[i * 3 + 2];
        Real cos_term = std::cos(2.0 * sctl::const_pi<Real>() * n * z / L);
        Real radicand = 1.0 + amplitude * cos_term * cos_omega_t;
        Real expected_radius = R0 * std::sqrt(radicand);
        Real expected_rdot = - (R0 * amplitude * omega * cos_term * sin_omega_t) / (2.0 * std::sqrt(radicand));
        Real actual_radius = wall.radiusIn()[i];
        Real actual_rdot = wall.rdotIn()[i];
        std::cout << std::setprecision(10);
        std::cout << "     Index " << i << " (z=" << z << "): "
                  << "Expected R=" << expected_radius << ", Got R=" << actual_radius
                  << "  |  Expected Rdot=" << expected_rdot << ", Got Rdot=" << actual_rdot
                  << std::endl;
        ASSERT_NEAR(actual_radius, expected_radius, 1e-10);
        ASSERT_NEAR(actual_rdot, expected_rdot, 1e-10);
    }
}
// ==========================================
// Test Suite & Main
// ==========================================

TEST_SUITE(func_wall_suite)
{
    TEST(constant_walls<double>);
    TEST(spatial_walls<double>);
    TEST(spatiotemporal_walls<double>);
    TEST(r2_spatiotemporal_walls<double>);
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