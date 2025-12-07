#ifndef __FUNC_WALL__
#include "../FuncWall.hpp"
#endif

template <class Real, class InnerFunc, class OuterFunc>
FuncWall<Real, InnerFunc, OuterFunc>::FuncWall( const Real dt, 
        const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
        const sctl::Vector<Real>& center_in,
        const sctl::Vector<Real>& radius_out,
        const sctl::Vector<Real>& radius_in,
        const InnerFunc& inner_func,
        const OuterFunc& outer_func)
    :   Wall<Real>(dt, center_out, center_in, radius_out, radius_in),
        _inner_func(inner_func),
        _outer_func(outer_func) {}


template <class Real, class InnerFunc, class OuterFunc>
FuncWall<Real, InnerFunc, OuterFunc>::~FuncWall() {}


template <class Real, class InnerFunc, class OuterFunc>
void FuncWall<Real, InnerFunc, OuterFunc>::update() {
    const Real t = this->getTime();

    // Inner wall
    if constexpr (std::is_invocable_v<
                      InnerFunc,
                      sctl::Vector<Real>&,      // radius
                      sctl::Vector<Real>&,      // rdot
                      Real,                     // time
                      const sctl::Vector<Real>& // coords
                  >)
    {
        // InnerFunc(radius, rdot, time, coords)
        _inner_func(this->_radius_in,
                    this->_rdot_in,
                    t,
                    this->_center_coords_in);
    }
    else if constexpr (std::is_invocable_v<
                           InnerFunc,
                           sctl::Vector<Real>&,
                           sctl::Vector<Real>&,
                           Real
                       >)
    {
        // InnerFunc(radius, rdot, time)
        _inner_func(this->_radius_in,
                    this->_rdot_in,
                    t);
    }
    else if constexpr (std::is_invocable_v<
                           InnerFunc,
                           sctl::Vector<Real>&,
                           Real,
                           const sctl::Vector<Real>&
                       >)
    {
        // InnerFunc(radius, time, coords) – no derivative
        _inner_func(this->_radius_in, t, this->_center_coords_in);
        this->_rdot_in = 0; // will this work?
    }
    else if constexpr (std::is_invocable_v<
                           InnerFunc,
                           sctl::Vector<Real>&,
                           Real
                       >)
    {
        // InnerFunc(radius, time) – no derivative
        _inner_func(this->_radius_in, t);
        this->_rdot_in = 0;
    }
    else {
        static_assert(sizeof(InnerFunc) == 0,
                      "InnerFunc has unsupported signature");
    }

    // =========================
    // Outer wall
    // =========================
    if constexpr (std::is_invocable_v<
                      OuterFunc,
                      sctl::Vector<Real>&,      // radius
                      sctl::Vector<Real>&,      // rdot
                      Real,                     // time
                      const sctl::Vector<Real>& // coords
                  >)
    {
        // OuterFunc(radius, rdot, time, coords)
        _outer_func(this->_radius_out,
                    this->_rdot_out,
                    t,
                    this->_center_coords_out);
    }
    else if constexpr (std::is_invocable_v<
                           OuterFunc,
                           sctl::Vector<Real>&,
                           sctl::Vector<Real>&,
                           Real
                       >)
    {
        // OuterFunc(radius, rdot, time)
        _outer_func(this->_radius_out,
                    this->_rdot_out,
                    t);
    }
    else if constexpr (std::is_invocable_v<
                           OuterFunc,
                           sctl::Vector<Real>&,
                           Real,
                           const sctl::Vector<Real>&
                       >)
    {
        // OuterFunc(radius, time, coords) – no derivative
        _outer_func(this->_radius_out, t, this->_center_coords_out);
        this->_rdot_out = 0; // will this work?
    }
    else if constexpr (std::is_invocable_v<
                           OuterFunc,
                           sctl::Vector<Real>&,
                           Real
                       >)
    {
        // OuterFunc(radius, time) – no derivative
        _outer_func(this->_radius_out, t);
        this->_rdot_out = 0;
    }
    else {
        static_assert(sizeof(OuterFunc) == 0,
                      "OuterFunc has unsupported signature");
    }

    // Make sure wall is physical after update
    this->enforceGapGeometry();
    // Increment time step
    this->incrementTimeStep();
}
