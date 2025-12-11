#ifndef __FUNC_WALL__
#define __FUNC_WALL__

#include "Wall.hpp"
#include <type_traits>

/**
 * @brief FuncWall class template for wall dynamics defined by user-provided functions, inherits from Wall base class.
 * 
 * templated functors for inner and outer wall dynamics allow varying signatures, dynamics and parameters
 * 
 * The functor signatures supported are:
 * 1. void func(sctl::Vector<Real>& radius, sctl::Vector<Real>& rdot, Real time, sctl::Vector<Real>& coords)
 * 2. void func(sctl::Vector<Real>& radius, sctl::Vector<Real>& rdot, Real time)
 * 3. void func(sctl::Vector<Real>& radius, Real time, sctl::Vector<Real>& coords)  // rdot set to zero
 * 4. void func(sctl::Vector<Real>& radius, Real time)  // rdot set to zero 
 * The functors are expected to modify the radius (and optionally its time derivative) based on the current time and spatial coordinates.
 * 
 * @tparam Real 
 * @tparam InnerFunc 
 * @tparam OuterFunc 
 */
template <class Real, class InnerFunc, class OuterFunc>
class FuncWall: public Wall<Real>{
    private:
        InnerFunc _inner_func;
        OuterFunc _outer_func;

        // Helper function to flexibly call wall functors with varying signatures
        template <class Func>
        void apply_wall_logic(Func& func, 
                            sctl::Vector<Real>& radius, 
                            sctl::Vector<Real>& rdot, 
                            Real t, 
                            sctl::Vector<Real>& coords,
                            sctl::Vector<Real>& cdot);

    public:
        FuncWall( const Real dt, 
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
            const sctl::Vector<Real>& radius_in,
            const sctl::Vector<Real>& radius_out,
            const InnerFunc& inner_func,
            const OuterFunc& outer_func
        );
        ~FuncWall() override;

        void update() override;

};

#include "detail/FuncWall.tpp"
#endif