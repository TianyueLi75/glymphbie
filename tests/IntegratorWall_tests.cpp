#include "GlymphBIE.hpp"
#include "unit_test_framework.hpp"
#include <cmath>
#include <utility> 
#include <iostream>
#include <iomanip> // For std::setprecision

// ==========================================
// Functor Definitions
// ==========================================


template <class Real, size_t hillExp>
struct HillFunc {
    // 
    Real _k_pow; 

    // Constructor: Computes K^N once and stores it.
    HillFunc(const Real K_half ) {
        _k_pow = 1.0;
        
        for (size_t i = 0; i < hillExp; ++i) {
            _k_pow *= K_half;
        }
    }

    // The Functor
    Real operator()(Real x) const {
        // 1. Calculate x^N
        Real x_pow = 1.0;
        
        // The compiler sees 'N' is constant and unrolls this loop.
        // It becomes: x_pow = x * x * x ... (using only multiplication)
        for (size_t i = 0; i < hillExp; ++i) {
            x_pow *= x;
        }

        // 2. Return Hill equation: x^N / (x^N + K^N)
        return x_pow / (x_pow + _k_pow);
    }
};



template <class Real>
struct InnerWallFunctor {
    // --- Parameters ---
    Real _R0;       // Base Radius
    Real _alpha_V;  // Vasodilation sensitivity
    Real _V0;       // Baseline V
    
    // Gaussian Amplitude Parameters
    Real _ac;       // Max amplitude
    Real _Vc;       // Center V for max amplitude
    Real _inv_2sigma2; 
    
    // Wave Parameters
    Real _omega;    // Frequency
    Real _k;        // Wavenumber

    InnerWallFunctor(Real R0, Real alpha_V, Real V0, 
                     Real ac, Real Vc, Real sigma_c, 
                     Real omega, Real k)
        : _R0(R0), _alpha_V(alpha_V), _V0(V0), 
          _ac(ac), _Vc(Vc), _omega(omega), _k(k) 
    {
        _inv_2sigma2 = 1.0 / (2.0 * sigma_c * sigma_c);
    }

    Real operator()(Real z, Real V, Real G, Real t) const {

        Real diff = V - _Vc;
        Real gaussian = std::exp(-(diff * diff) * _inv_2sigma2);
        Real Ac = _ac * gaussian;


        Real s = z; 
        Real phase = (_omega * t) + (_k * s);
        Real sin_val = std::sin(phase);


        Real term_linear = 1.0 + _alpha_V * (V - _V0);
        Real term_wave = Ac * sin_val;

        return _R0 * (term_linear + term_wave);
    }
};

template <class Real>
struct OuterWallFunctor {
    Real _R0;
    Real _alpha_G;
    Real _G0;

    OuterWallFunctor(Real R0, Real alpha_G, Real G0)
        : _R0(R0), _alpha_G(alpha_G), _G0(G0) {}

    Real operator()(Real z, Real V, Real G, Real t) const {
        Real term_swell = 1.0 - _alpha_G * (G - _G0);
        return _R0 * term_swell;
    }
};

// ==========================================
// Test Cases
// ==========================================
template <class Real>
TEST_CASE(constant_walls){
    using HillV_t  = HillFunc<Real, 2>;
    using HillG_t  = HillFunc<Real, 2>;
    using Inner_t  = InnerWallFunctor<Real>;
    using Outer_t  = OuterWallFunctor<Real>;
    using Solver_t = IntegratorWall<Real, HillV_t, HillG_t, Inner_t, Outer_t>;

    // set parameters
    Real dt = 0.1;
    Real R_in_base = 1.0;
    Real R_out_base = 1.2;

    Real alpha_V = 0.0;
    Real V0 = 0.0; 
    Real alpha_G = 0.0;
    Real G0 = 0.0;

    Real ac = 0.0; Real Vc = 0.5; Real sigma = 0.1; Real omega = 0; Real k = 0;

    // instantiate functors
    HillV_t f_V(0.5);
    HillG_t f_G(0.5);
    Inner_t f_inner(R_in_base, alpha_V, V0, ac, Vc, sigma, omega, k);
    Outer_t f_outer(R_out_base, alpha_G, G0);
    // setup data vectors
    long N = 10;
    sctl::Vector<Real> c_in(N*3), c_out(N*3);
    sctl::Vector<Real> r_in(N), r_out(N);
    sctl::Vector<Real> activity(N), V(N), G(N);

    // 
    c_in = 0; c_out = 0;
    for (long i=0; i<N; ++i) {
        r_in[i]  = R_in_base;
        r_out[i] = R_out_base;
        activity[i] = 1.0; 
        V[i]   = V0; 
        G[i]   = G0; 
    }
    // instantiate integrator
    Real tau = 1.0; 
    Solver_t solver(
        dt, c_out, c_in, r_out, r_in, 
        activity, V, G, tau, tau, 
        f_V, f_G, f_inner, f_outer
    ); 
    solver.update();
    // verify radii remain constant
    for (long i=0; i<N; ++i) {
        ASSERT_NEAR(solver.radiusIn()[i], R_in_base, 1e-10);
        ASSERT_NEAR(solver.radiusOut()[i], R_out_base, 1e-10);
    }


}


// template <class Real>
// TEST_CASE(bio_parameters) {
//     using HillV_t  = HillFunc<Real, 2>;
//     using HillG_t  = HillFunc<Real, 2>;
//     using Inner_t  = InnerWallFunctor<Real>;
//     using Outer_t  = OuterWallFunctor<Real>;
//     using Solver_t = IntegratorWall<Real, HillV_t, HillG_t, Inner_t, Outer_t>;

//     // set parameters
//     Real dt = 0.1;
//     Real R_in_base = 1.0;
//     Real R_out_base = 1.2;
// }


// =========================================================
// The Generic Test Case
// =========================================================
template <class Real>
TEST_CASE(inner_outer_wall_integration) {
    
    // 1. Type Definitions
    // -----------------------------------------------------
    using HillV_t  = HillFunc<Real, 4>; // Deg 4
    using HillG_t  = HillFunc<Real, 2>; // Deg 2
    using Inner_t  = InnerWallFunctor<Real>;
    using Outer_t  = OuterWallFunctor<Real>;
    using Solver_t = IntegratorWall<Real, HillV_t, HillG_t, Inner_t, Outer_t>;

    // 2. Physical Parameters
    // -----------------------------------------------------
    Real dt = 0.1;
    
    // Functor Params
    Real R_in_base = 1.0;
    Real R_out_base = 1.15; // Starts with 0.15 gap
    
    // Inner Wall: High sensitivity to V to force rapid expansion
    Real alpha_V = 1.0; 
    Real V0 = 0.0;
    // Turn off waves/gaussian for simplicity to isolate expansion logic
    Real ac = 0.0; Real Vc = 0.5; Real sigma = 0.1; Real omega = 0; Real k = 0;

    // Outer Wall: Static (alpha_G = 0)
    // It should NOT move unless forced by the gap constraint
    Real alpha_G = 0.0; 
    Real G0 = 0.0;

    // 3. Instantiate Functors
    // -----------------------------------------------------
    HillV_t f_V(0.5);
    HillG_t f_G(0.5);
    Inner_t f_inner(R_in_base, alpha_V, V0, ac, Vc, sigma, omega, k);
    Outer_t f_outer(R_out_base, alpha_G, G0);

    // 4. Setup Data Vectors (1 Particle)
    // -----------------------------------------------------
    long N = 1;
    sctl::Vector<Real> c_in(N*3), c_out(N*3);
    sctl::Vector<Real> r_in(N), r_out(N);
    sctl::Vector<Real> act(N), V(N), G(N);

    // Initial Conditions
    c_in = 0; c_out = 0; // Centers aligned (delta = 0)
    
    // Set Center Z to 0.0 for the functor lookup
    c_in[2] = 0.0; c_out[2] = 0.0;

    r_in[0]  = R_in_base;  // 1.0
    r_out[0] = R_out_base; // 1.15

    act[0] = 1.0; 
    V[0]   = 0.0; 
    G[0]   = 0.0; 

    // 5. Instantiate Integrator
    // -----------------------------------------------------
    Real tau = 1.0; 
    
    Solver_t solver(
        dt, c_out, c_in, r_out, r_in, 
        act, V, G, tau, tau, 
        f_V, f_G, f_inner, f_outer
    );
    

    solver.update();

    // 7. Verify Logic
    // -----------------------------------------------------
    
    // Retrieve new values
    Real r_in_new  = solver.radiusIn()[0];
    Real r_out_new = solver.radiusOut()[0];
    Real v_in_new  = solver.rdotIn()[0];
    Real v_out_new = solver.rdotOut()[0];

    // EXPECTATION 1: Inner wall expanded
    // V increased -> Radius = R0 * (1 + alpha*V) -> Expanded
    SCTL_ASSERT(r_in_new > R_in_base);

    // EXPECTATION 2: Finite Difference Velocity (Inner)
    // v = (r_new - r_old) / dt
    Real expected_v_in = (r_in_new - R_in_base) / dt;
    ASSERT_NEAR(v_in_new, expected_v_in, 1e-5);

    // EXPECTATION 3: Gap Constraint & Outer Velocity
    Real default_min_gap = 1e-5; // Based on your header defaults

    if (r_in_new + default_min_gap > R_out_base) {
        std::cout << "Constraint Activated (Good)." << std::endl;

        // Check Position: r_out must be exactly r_in + gap
        Real expected_r_out = r_in_new + default_min_gap;
        ASSERT_NEAR(r_out_new, expected_r_out, 1e-5);

        // Check Velocity: Outer wall moved due to collision
        // It should have positive velocity even though alpha_G = 0
        Real expected_v_out = (r_out_new - R_out_base) / dt;
        ASSERT_NEAR(v_out_new, expected_v_out, 1e-5);

    } else {
        std::cout << "Constraint Not Activated (Inner didn't expand enough)." << std::endl;
        // If constraint didn't activate, outer wall should be static
        ASSERT_NEAR(v_out_new, 0.0, 1e-5);
    }
}

TEST_SUITE(Integrator_tests)
{
    using Real = double;
    TEST(inner_outer_wall_integration<Real>);
    TEST(constant_walls<Real>);
}


int main(int argc, char** argv) {
    std::cout << "=========================================" << std::endl;
    std::cout << " STARTING INTEGRATORWALL UNIT TESTS" << std::endl;
    std::cout << "=========================================" << std::endl;
    
    RUN_SUITE(Integrator_tests);
    
    std::cout << "\n=========================================" << std::endl;
    std::cout << " ALL TESTS COMPLETED" << std::endl;
    std::cout << "=========================================" << std::endl;
    return 0;
}

