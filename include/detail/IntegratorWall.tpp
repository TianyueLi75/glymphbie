#ifndef __INTEGRATOR_WALL__
#include "../IntegratorWall.hpp"
#endif
template <class Real, class NonLinearFuncV, class NonLinearFuncG, class InnerFunc, class OuterFunc>
IntegratorWall<Real, NonLinearFuncV, NonLinearFuncG, InnerFunc, OuterFunc>::IntegratorWall( const Real dt, 
        const sctl::Vector<Real>& center_in,
        const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
        const sctl::Vector<Real>& radius_in,
        const sctl::Vector<Real>& radius_out,
        const sctl::Vector<Real>& neuronal_activity,
        const sctl::Vector<Real>& V,
        const sctl::Vector<Real>& G,
        const Real tau_V,
        const Real tau_G,
        const NonLinearFuncV& FV,
        const NonLinearFuncG& FG,
        const InnerFunc& inner_func,
        const OuterFunc& outer_func
        )

    :   Wall<Real>(dt, center_in, center_out, radius_in, radius_out),
        _neuronal_activity(neuronal_activity),
        _V(V),
        _G(G),
        _tau_V(tau_V),
        _tau_G(tau_G),
        _FV(FV),
        _FG(FG),
        _inner_func(inner_func),
        _outer_func(outer_func) {}


template <class Real, class NonLinearFuncV, class NonLinearFuncG, class InnerFunc, class OuterFunc>
IntegratorWall<Real, NonLinearFuncV, NonLinearFuncG, InnerFunc, OuterFunc>::~IntegratorWall() {}


template <class Real, class NonLinearFuncV, class NonLinearFuncG, class InnerFunc, class OuterFunc>
void IntegratorWall<Real, NonLinearFuncV, NonLinearFuncG, InnerFunc, OuterFunc>::expIntStep(){
    const Real dt = this->_dt;
    const Real dt_over_one = 1.0 / dt;
    const Real t_next = this->getTime() + dt;
    const Real decayV = std::exp(-dt/_tau_V);
    const Real decayG = std::exp(-dt/_tau_G);

    for (size_t i = 0; i < this->_V.Dim(); i++) {
        const Real N = _neuronal_activity[i];
        const Real V_inf = _FV(N);
        const Real G_inf = _FG(N);
        const Real offsetY = this->_center_coords_in[3*i+1] - this->_center_coords_out[3*i+1];
        const Real offsetZ = this->_center_coords_in[3*i+2] - this->_center_coords_out[3*i+2];
        const Real delta = std::sqrt(offsetZ*offsetZ + offsetY*offsetY);
        
        // exponential integration step
        _V[i] = V_inf + (_V[i] - V_inf)*decayV;
        _G[i] = G_inf + (_G[i] - G_inf)*decayG;

        // capture old radius values for velocity calculation
        Real rin_old = this->_radius_in[i];
        Real rout_old = this->_radius_out[i];

        // compute new physical radii based on wall functions
        Real rin_phys = _inner_func(this->_center_coords_in[3*i], _V[i], _G[i], t_next);
        Real rout_phys = _outer_func(this->_center_coords_out[3*i], _V[i], _G[i], t_next);

        assert(rin_phys > 0.0 && rout_phys > 0.0);
        const Real minimum_gap = rout_phys - rin_phys - delta;
        if (minimum_gap < this->_min_gap) {
            // Adjust outer radius to enforce minimum gap
            rout_phys = rin_phys + delta + this->_min_gap;
        }

        Real rin_dot = (rin_phys - rin_old) * dt_over_one;
        Real rout_dot = (rout_phys - rout_old) * dt_over_one;

        // Update wall radii and velocities
        this->_radius_in[i] = rin_phys;
        this->_radius_out[i] = rout_phys;
        this->_rdot_in[i] = rin_dot;
        this->_rdot_out[i] = rout_dot;

    }
    

}

template <class Real, class NonLinearFuncV, class NonLinearFuncG, class InnerFunc, class OuterFunc>
void IntegratorWall<Real, NonLinearFuncV, NonLinearFuncG, InnerFunc, OuterFunc>::update() {
    // Perform exponential integration step to update wall radii and velocities
    this->expIntStep();
    // Increment time step
    this->incrementTimeStep();

}
