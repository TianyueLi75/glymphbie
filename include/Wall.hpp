#ifndef __GLYMPH_WALL__
#define __GLYMPH_WALL__

#include <csbq.hpp>
#include <stdexcept>

/**
 * @brief Wall class template for wall dynamics in glymphatic flow simulations.
 * 
 * @tparam Real 
 * @param dt Time step size for wall dynamics.
 * @param center_out Vector of outer wall center coordinates (size: 3 * number of grid points).
 * @param center_in Vector of inner wall center coordinates (size: 3 * number of grid points).
 * @param radius_out Vector of outer wall radii (size: number of grid points).
 * @param radius_in Vector of inner wall radii (size: number of grid points).
 * @note The wall dynamics (update of radii and their time derivatives) are to be defined in derived classes.
 */
template <class Real>
class Wall
{
    protected:
        Real _dt;
        Real _min_gap = 1e-5;
        Real _min_inner_radius = 1e-5;
        size_t _time_step = 0;
        sctl::Vector<Real> _center_coords_out;
        sctl::Vector<Real> _center_coords_in;
        sctl::Vector<Real> _cdot_in;
        sctl::Vector<Real> _cdot_out; 
        sctl::Vector<Real> _radius_out;
        sctl::Vector<Real> _radius_in;
        sctl::Vector<Real> _rdot_out;
        sctl::Vector<Real> _rdot_in;
        size_t _spatial_grid_size = 0;

    public:
        Wall(const Real dt, 
            const sctl::Vector<Real>& center_out, // Spatial grid for inner and outer walls is equivalent to the x coordinates of the centerlines
            const sctl::Vector<Real>& center_in,
            const sctl::Vector<Real>& radius_out,
            const sctl::Vector<Real>& radius_in
            // add sctl::Comm if needed for MPI
        );
        virtual ~Wall();

        virtual void update() = 0;

        void enforceGapGeometry();
        
        // Replace grid and data in case mesh is refined
        virtual void setGridAndData(
                                const sctl::Vector<Real>& center_out,
                                const sctl::Vector<Real>& center_in,
                                const sctl::Vector<Real>& radius_out,
                                const sctl::Vector<Real>& radius_in, 
                                const sctl::Vector<Real>& rdot_out,
                                const sctl::Vector<Real>& rdot_in
                                );

        



        // Accessors
        void incrementTimeStep() { _time_step++; }
        size_t getTimeStep() const { return _time_step; }
        Real getTime() const { return _time_step*_dt; }
        const sctl::Vector<Real>& centerCoordsOut() const { return _center_coords_out; }
        const sctl::Vector<Real>& centerCoordsIn() const { return _center_coords_in; }
        const sctl::Vector<Real>& cdotOut() const { return _cdot_out; }
        const sctl::Vector<Real>& cdotIn() const { return _cdot_in; }
        const sctl::Vector<Real>& radiusOut() const { return _radius_out; }
        const sctl::Vector<Real>& radiusIn() const { return _radius_in; }
        const sctl::Vector<Real>& rdotOut() const { return _rdot_out; }
        const sctl::Vector<Real>& rdotIn() const { return _rdot_in; }
};
#include "detail/Wall.tpp"
#endif // __GLYMPH_WALL__