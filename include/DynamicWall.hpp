#ifndef __DYNAMIC_WALL__
#define __DYNAMIC_WALL__

#include "Wall.hpp"
#include <type_traits>
#include <cmath>
#include <cassert>


template <class Real, class LinearFunc, class NonLinearFunc, class InnerFunc, class OuterFunc>
class DynamicWall: public Wall<Real>{
    private:
        sctl::Vector<Real> _Y; // state vector
        LinearFunc _A; // linear operator
        NonLinearFunc _F; // non-linear operator
        InnerFunc _inner_func;
        OuterFunc _outer_func;

    public:
        DynamicWall( const Real dt, 
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& center_out, 
            const sctl::Vector<Real>& radius_in,
            const sctl::Vector<Real>& radius_out,
            const sctl::Vector<Real>& initial_state,
            const LinearFunc& A,
            const NonLinearFunc& F,
            const InnerFunc& inner_func,
            const OuterFunc& outer_func
        );
        ~DynamicWall() override;

        void update() override;

        void expIntStep();

        

};

#include "detail/DynamicWall.tpp"
#endif