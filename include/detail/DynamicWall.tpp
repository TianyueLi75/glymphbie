#ifndef __DYNAMIC_WALL__
#include "../DynamicWall.hpp"
#endif
template <class Real, class LinearFunc, class NonLinearFunc, class InnerFunc, class OuterFunc>
DynamicWall<Real, LinearFunc, NonLinearFunc,  InnerFunc, OuterFunc>::DynamicWall( const Real dt, 
        const sctl::Vector<Real>& center_in,
        const sctl::Vector<Real>& center_out, 
        const sctl::Vector<Real>& radius_in,
        const sctl::Vector<Real>& radius_out,
        const sctl::Vector<Real>& Y0,
        const LinearFunc& A,
        const NonLinearFunc& F,
        const InnerFunc& inner_func,
        const OuterFunc& outer_func
        )

    :   Wall<Real>(dt, center_in, center_out, radius_in, radius_out),
        _Y(Y0),
        _A(A),
        _F(F),
        _inner_func(inner_func),
        _outer_func(outer_func)    
    {}


template <class Real, class LinearFunc, class NonLinearFunc, class InnerFunc, class OuterFunc>
DynamicWall<Real, LinearFunc, NonLinearFunc,  InnerFunc, OuterFunc>::~DynamicWall() {}


template <class Real, class LinearFunc, class NonLinearFunc, class InnerFunc, class OuterFunc>
void DynamicWall<Real, LinearFunc, NonLinearFunc,  InnerFunc, OuterFunc>::expIntStep(){
    const Real dt = this->_dt;
    

    for (size_t i = 0; i < this->_spatial_grid_size; i++) {
    
    }
    

}

template <class Real, class LinearFunc, class NonLinearFunc, class InnerFunc, class OuterFunc>
void DynamicWall<Real, LinearFunc, NonLinearFunc,  InnerFunc, OuterFunc>::update() {
    // Perform exponential integration step to update wall radii and velocities
    this->expIntStep();
    // Increment time step
    this->incrementTimeStep();

}
