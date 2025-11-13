#ifndef __FUNC_WALL__
#define __FUNC_WALL__

#include "Wall.hpp"
#include <functional>

template <class Real>
class FuncWall: public Wall<Real>{
    private:
        std::function<void(sctl::Vector<Real>&, Real)> _inner_func;
        std::function<void(sctl::Vector<Real>&, Real)> _outer_func;

    public:
        FuncWall(size_t n_tsteps, Real dt, size_t grid_size, 
        std::function<void(sctl::Vector<Real>&, Real)> inner_func,
        std::function<void(sctl::Vector<Real>&, Real)> outer_func
        );
        ~FuncWall() override;

        void update() override;

};
#endif