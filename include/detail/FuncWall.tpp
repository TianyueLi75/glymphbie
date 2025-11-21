
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
    Real t = this->getTime();

    if constexpr (std::is_invocable_v<InnerFunc, sctl::Vector<Real>&, Real, const sctl::Vector<Real>&>) {
        // InnerFunc is callable with (sctl::Vector<Real>&, Real , const sctl::Vector<Real>&)
        _inner_func(this->_radius_in, t, this->_center_coords_in);
    } else {
        // Assume InnerFunc is callable with (sctl::Vector<Real>&, Real)
        _inner_func(this->_radius_in, t);
    }

    if constexpr (std::is_invocable_v<OuterFunc, sctl::Vector<Real>&, Real, const sctl::Vector<Real>&>) {
        // OuterFunc is callable with (sctl::Vector<Real>&, Real , const sctl::Vector<Real>&)
        _outer_func(this->_radius_out, t, this->_center_coords_out);
    } else {
        // Assume OuterFunc is callable with (sctl::Vector<Real>&, Real)
        _outer_func(this->_radius_out, t);
    }

    // Make sure wall is physical after update
    this->enforceGapGeometry(); 
    // Increment time step
    this->incrementTimeStep();

}
