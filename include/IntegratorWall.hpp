#ifndef __INTEGRATOR_WALL__
#define __INTEGRATOR_WALL__

#include "Wall.hpp"
#include <type_traits>
#include <cmath>
#include <cassert>
// IntegratorWall class template work for wall dynamics defined by vasodilation and glial swelling variables
// with non-linear functions defining their steady-state values based on neuronal activity input which is time/space dependent, but does not depend on G,V
/**
 * @brief IntegratorWall class template for wall dynamics defined by vasodilation and glial swelling variables, inherits from Wall base class.
 * 
 * The wall radii are updated based on internal state variables V (vasodilation) and G (glial swelling) which evolve according to exponential integration towards steady-state values defined by non-linear functions of neuronal activity.
 * The non-linear functions FV and FG define the steady-state values of V and G respectively based on the neuronal activity input.
 * The wall dynamics are further defined by user-provided functors for inner and outer wall radius calculations based on the current V, G, time, and spatial coordinates.
 * 
 * @tparam Real 
 * @tparam NonLinearFuncV 
 * @tparam NonLinearFuncG 
 * @tparam InnerFunc 
 * @tparam OuterFunc 
 */
template <class Real, class NonLinearFuncV, class NonLinearFuncG, class InnerFunc, class OuterFunc>
class IntegratorWall: public Wall<Real>{
    private:
        sctl::Vector<Real> _neuronal_activity;
        sctl::Vector<Real> _V; // vasodilation variable [0,1]
        sctl::Vector<Real> _G; // glial swelling variable [0,1]
        Real _tau_V; // time constant for vasodilation
        Real _tau_G; // time constant for glial swelling
        NonLinearFuncV _FV; // Hill function for vasodilation
        NonLinearFuncG _FG; // Hill function for glial swelling  
        InnerFunc _inner_func;
        OuterFunc _outer_func;

    public:
        IntegratorWall( const Real dt, 
            const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& radius_out,
            const sctl::Vector<Real>& radius_in,
            const sctl::Vector<Real>& neuronal_activity,
            const sctl::Vector<Real>& V,
            const sctl::Vector<Real>& G,
            const Real tau_V,
            const Real tau_G,
            const NonLinearFuncV& FV,
            const NonLinearFuncG& FG,
            const InnerFunc& inner_func,
            const OuterFunc& outer_func
        );
        ~IntegratorWall() override;

        void update() override;

        void expIntStep();

        const sctl::Vector<Real>& getNeuronalActivity() const { return _neuronal_activity; }
        

};

#include "detail/IntegratorWall.tpp"
#endif