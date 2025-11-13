#include "Wall.hpp"
#include <stdexcept>

template <class Real>
Wall<Real>::Wall(const size_t n_tsteps, const Real dt, const size_t grid_size)
    : _n_tsteps(n_tsteps), _dt(dt), _grid_size(grid_size),
      _center_coords_out(3 * grid_size),
      _center_coords_in(3 * grid_size),
      _radius_out(grid_size),
      _radius_in(grid_size) {}

template <class Real>
Wall<Real>::~Wall() {}

template <class Real>
void Wall<Real>::setGridAndData(size_t new_grid_size,
                                   const sctl::Vector<Real>& center_out,
                                   const sctl::Vector<Real>& center_in,
                                   const sctl::Vector<Real>& radius_out,
                                   const sctl::Vector<Real>& radius_in) {
    if (center_out.Dim() != 3 * new_grid_size ||
        center_in.Dim() != 3 * new_grid_size ||
        radius_out.Dim() != new_grid_size ||
        radius_in.Dim() != new_grid_size) {
        throw std::invalid_argument("Input sizes do not match new grid size");
    }
    _grid_size = new_grid_size;
    _center_coords_out = center_out;
    _center_coords_in = center_in;
    _radius_out = radius_out;
    _radius_in = radius_in;
}

// Explicit instantiation
template class Wall<float>;
template class Wall<double>;