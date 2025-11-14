#include "Wall.hpp"
#include <stdexcept>
#include <cmath>
#include <cassert>
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

template <class Real>
void Wall<Real>::enforceGapGeometry() {
    // TODO: add better checking for negative radius
    for (size_t i = 0; i < _grid_size; ++i){
        const size_t j = 3*i;

        const Real deltaX = _center_coords_out[j] - _center_coords_in[j];
        const Real deltaY = _center_coords_out[j+1] - _center_coords_in[j+1];
        const Real deltaZ = _center_coords_out[j+2] - _center_coords_in[j+2];
        // deltaX should be zero for planar walls, but just in case. Its worth considering how things may change for non-planar walls.
        const Real dist = std::sqrt(deltaX*deltaX + deltaY*deltaY+ deltaZ*deltaZ);

        Real& r_out = _radius_out[i];
        Real& r_in = _radius_in[i];

        // Enforce positive radii and outer radius larger than inner radius
        assert(r_in > 0 && r_out > r_in);

        // Enforce minimum inner radius
        if (r_in < _min_inner_radius){
            r_in = _min_inner_radius;
            std::cerr << "Warning: inner radius below minimum at x = " << _center_coords_in[j] << ", adjusting inner radius to enforce minimum." << std::endl;
        }
        // Enforce minimum gap by bumping outer radius
        const Real gap = r_out - r_in - dist;
        if (gap < _min_gap){
            r_out = dist + r_in + _min_gap;
            std::cerr << "Warning: gap below minimum at x = " << _center_coords_out[j] << ", adjusting outer radius to enforce minimum gap." << std::endl;
        }
    }
}
// Explicit instantiation
template class Wall<float>;
template class Wall<double>;