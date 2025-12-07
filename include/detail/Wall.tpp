#ifndef __WALL__
#include "../Wall.hpp"
#endif

#include <cmath>
#include <iostream> // For std::cerr, std::endl
#include <cassert>  // For assert()


template <class Real>
Wall<Real>::Wall(const Real dt, 
            const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& radius_out,
            const sctl::Vector<Real>& radius_in)
    : _dt(dt), _center_coords_out(center_out), _center_coords_in(center_in),
      _radius_out(radius_out), _radius_in(radius_in),
      _rdot_in(radius_in.Dim()),
      _rdot_out(radius_out.Dim())

      
{
    // Validation checks
    if (center_out.Dim() != radius_out.Dim() * 3 || center_in.Dim() != radius_in.Dim() * 3) {
        throw std::invalid_argument("Center coordinate vector size must be 3 times the radius vector size.");
    }
    // for now, require same spatial grid for inner and outer walls if not the enforce radius function won't work
    // note this doesn't check actual coordinates, just vector sizes
    if (center_out.Dim() != center_in.Dim()) {
        throw std::invalid_argument("Outer and inner center coordinate vectors must have the same size.");
    }
    if (radius_out.Dim() != radius_in.Dim()) {
        throw std::invalid_argument("Outer and inner radius vectors must have the same size.");
    }
    _spatial_grid_size = radius_out.Dim();

    _rdot_in.SetZero();
    _rdot_out.SetZero();
}

template <class Real>
Wall<Real>::~Wall() {}

template <class Real>
void Wall<Real>::setGridAndData(
                                   const sctl::Vector<Real>& center_out,
                                   const sctl::Vector<Real>& center_in,
                                   const sctl::Vector<Real>& radius_out,
                                   const sctl::Vector<Real>& radius_in,
                                   const sctl::Vector<Real>& rdot_out,
                                   const sctl::Vector<Real>& rdot_in
                                   ) {
    if (center_out.Dim() != radius_out.Dim()*3 || center_in.Dim() != radius_in.Dim()*3) {
        throw std::invalid_argument("in/out radius and in/out center coordinate vectors must match grid size.");
    }
    if (center_out.Dim() != center_in.Dim()) {
        throw std::invalid_argument("inner and outer center coordinate vectors must have the same size for now");
    }
    _spatial_grid_size = radius_out.Dim();
    _center_coords_out = center_out;
    _center_coords_in = center_in;
    _radius_out = radius_out;
    _radius_in = radius_in;
    _rdot_out = rdot_out;
    _rdot_in = rdot_in;
}

template <class Real>
void Wall<Real>::enforceGapGeometry() {
    // TODO: add better checking for negative radius
    // TODO: what if spatial grid is different for inner and outer walls?
    for (size_t i = 0; i < _spatial_grid_size; ++i){
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