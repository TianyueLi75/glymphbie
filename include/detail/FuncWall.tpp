#ifndef __FUNC_WALL__
#include "../FuncWall.hpp"
#endif

template <class Real, class InnerFunc, class OuterFunc>
FuncWall<Real, InnerFunc, OuterFunc>::FuncWall( const Real dt, 
        const sctl::Vector<Real>& center_in,
        const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
        const sctl::Vector<Real>& radius_in,
        const sctl::Vector<Real>& radius_out,
        const InnerFunc& inner_func,
        const OuterFunc& outer_func)
    :   Wall<Real>(dt, center_in, center_out, radius_in, radius_out),
        _inner_func(inner_func),
        _outer_func(outer_func) {}


template <class Real, class InnerFunc, class OuterFunc>
FuncWall<Real, InnerFunc, OuterFunc>::~FuncWall() {}

template <class Real, class InnerFunc, class OuterFunc>

// adapt to approximate derivative with finite differences if needed
template <class Func> 
void FuncWall<Real, InnerFunc, OuterFunc>::apply_wall_logic(
    Func& func, 
    sctl::Vector<Real>& radius, 
    sctl::Vector<Real>& rdot, 
    Real t, 
    sctl::Vector<Real>& coords,
    sctl::Vector<Real>& cdot) 
{   
    // Check: (radius, rdot, coords, cdot, time) - covers both const& and & coords
    if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, sctl::Vector<Real>&, sctl::Vector<Real>&, sctl::Vector<Real>&, Real>) {
        func(radius, rdot, coords, cdot, t);
    }

    // Check: (radius, rdot, coords, time) - covers both const& and & coords
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, sctl::Vector<Real>&, sctl::Vector<Real>&, Real>) {
        func(radius, rdot, coords, t);
        cdot = 0.;
    }
    // Check: (radius, rdot, coords, time)
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, sctl::Vector<Real>&, Real>) {
        func(radius, rdot, coords, t);
        cdot = 0.;
    }
    // Check: (radius, coords, time) - No derivative
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, sctl::Vector<Real>&, Real>) {
        func(radius, coords, t);
        rdot = 0.;
        cdot = 0.;
    }
    // Check: (radius, coords) - No derivative
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, sctl::Vector<Real>&>) {
        func(radius, coords);
        rdot = 0.;
        cdot = 0.;
    }
    // Check: (radius, time) - No derivative
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&, Real>) {
        func(radius, t);
        rdot = 0.;
        cdot = 0.;
    }
    else if constexpr (std::is_invocable_v<Func, sctl::Vector<Real>&>) {
        func(radius);
        rdot = 0.;
        cdot = 0.;
    }
    else {
        static_assert(sizeof(Func) == 0, "Wall functor has unsupported signature");
    }
}





// need to be careful about when we evaluate fluid/particles vs wall update to ensure consistency right now wall functional form + implied volume conservation assumes update happens before fluid/particle eval at each time step 
template <class Real, class InnerFunc, class OuterFunc>
void FuncWall<Real, InnerFunc, OuterFunc>::update() {
    const Real t = this->getTime();

    // 1. Update Inner Wall
    apply_wall_logic(_inner_func, 
                     this->_radius_in, 
                     this->_rdot_in, 
                     t, 
                     this->_center_coords_in, this->_cdot_in);

    // 2. Update Outer Wall
    apply_wall_logic(_outer_func, 
                     this->_radius_out, 
                     this->_rdot_out, 
                     t, 
                     this->_center_coords_out, this->_cdot_out);

    // Make sure wall is physical after update
    this->enforceGapGeometry();
    // Increment time step
    this->incrementTimeStep();
}

template <class Real, class InnerFunc, class OuterFunc>
void FuncWall<Real, InnerFunc, OuterFunc>::initialize() {
    const Real t = this->getTime();
    // 1. Update Inner Wall
    apply_wall_logic(_inner_func, 
                     this->_radius_in, 
                     this->_rdot_in, 
                     t, 
                     this->_center_coords_in, this->_cdot_in);

    // 2. Update Outer Wall
    apply_wall_logic(_outer_func, 
                     this->_radius_out, 
                     this->_rdot_out, 
                     t, 
                     this->_center_coords_out, this->_cdot_out);

    // Make sure wall is physical after update
    this->enforceGapGeometry();
}
