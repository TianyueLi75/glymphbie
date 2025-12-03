
template <class Real, class HillV, class HillG, class InnerFunc, class OuterFunc>
IntegratorWall<Real, HillV, HillG, InnerFunc, OuterFunc>::IntegratorWall( const Real dt, 
        const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
        const sctl::Vector<Real>& center_in,
        const sctl::Vector<Real>& radius_out,
        const sctl::Vector<Real>& radius_in,
        const sctl::Vector<Real>& neuronal_activity,
        const sctl::Vector<Real>& V,
        const sctl::Vector<Real>& G,
        const Real tau_V,
        const Real tau_G,
        const HillV& hill_V,
        const HillG& hill_G,
        const InnerFunc& inner_func,
        const OuterFunc& outer_func
        )
    :   Wall<Real>(dt, center_out, center_in, radius_out, radius_in),
        _neuronal_activity(neuronal_activity),
        _V(V),
        _G(G),
        _tau_V(tau_V),
        _tau_G(tau_G),
        _hill_V(hill_V),
        _hill_G(hill_G),
        _inner_func(inner_func),
        _outer_func(outer_func) {}


template <class Real, class HillV, class HillG, class InnerFunc, class OuterFunc>
IntegratorWall<Real, HillV, HillG, InnerFunc, OuterFunc>::~IntegratorWall() {}


template <class Real, class HillV, class HillG, class InnerFunc, class OuterFunc>
void IntegratorWall<Real, HillV, HillG, InnerFunc, OuterFunc>::implicitStep(){
    this->getTimeStep();
    Real decay_V = this->_dt / this->_tau_V;
    Real decay_G = this->_dt / this->_tau_G;
    Real one_plus_decay_V = 1.0 + decay_V;
    Real one_plus_decay_G = 1.0 + decay_G;

    for (size_t i = 0; i < this->_spatial_grid_size; i++){

    }

}

template <class Real, class HillV, class HillG, class InnerFunc, class OuterFunc>
void IntegratorWall<Real, HillV, HillG, InnerFunc, OuterFunc>::update() {
    Real t = this->getTime();

    

    // Make sure wall is physical after update
    this->enforceGapGeometry(); 
    // Increment time step
    this->incrementTimeStep();

}
