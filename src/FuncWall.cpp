#include "FuncWall.hpp"
#include <stdexcept>
#include <iostream>

template <class Real>
FuncWall<Real>::FuncWall(size_t n_tsteps, Real dt, size_t grid_size,
        std::function<void(sctl::Vector<Real>&, Real)> inner_func,
        std::function<void(sctl::Vector<Real>&, Real)> outer_func)
    :   Wall<Real>(n_tsteps, dt, grid_size),
        _inner_func(inner_func),
        _outer_func(outer_func) {}


template <class Real>
FuncWall<Real>::~FuncWall() {}



template <class Real>
void FuncWall<Real>::update() {
    static size_t step = 0;
    // TODO: modify if we want it to be able to keep going
    if (step >= this->_n_tsteps) {
        std::cerr << "Warning: update called beyond time steps" << std::endl;
        return;
    }
    Real t = step * this->_dt;
    _inner_func(this->_radius_in, t);
    _outer_func(this->_radius_out, t);
    // TODO make sure wall is physical not just radii need to check w.r.t. center coordinates.
    step++;

}


// Explicit instantiation
template class FuncWall<float>;
template class FuncWall<double>;