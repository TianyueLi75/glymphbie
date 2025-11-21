#ifndef __FUNC_WALL__
#define __FUNC_WALL__

#include "Wall.hpp"
#include <type_traits>

template <class Real, class InnerFunc, class OuterFunc>
class FuncWall: public Wall<Real>{
    private:
        InnerFunc _inner_func;
        OuterFunc _outer_func;

    public:
        FuncWall( const Real dt, 
            const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& radius_out,
            const sctl::Vector<Real>& radius_in,
            const InnerFunc& inner_func,
            const OuterFunc& outer_func
        );
        ~FuncWall() override;

        void update() override;

};

#include "detail/FuncWall.tpp"
#endif