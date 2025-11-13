#ifndef __GLYMPH_WALL__
#define __GLYMPH_WALL__

#include <sctl.hpp>




template <class Real>
class Wall
{
    protected:
        size_t _n_tsteps;
        Real _dt;
        size_t _grid_size;

        sctl::Vector<Real> _center_coords_out;
        sctl::Vector<Real> _center_coords_in;

        sctl::Vector<Real> _radius_out;
        sctl::Vector<Real> _radius_in;

    public:
        Wall(const size_t n_tsteps, const Real dt, const size_t grid_size);
        virtual ~Wall();

        virtual void update() = 0;


        // Replace grid and data in case mesh is refined
        virtual void setGridAndData(size_t new_grid_size,
                                const sctl::Vector<Real>& center_out,
                                const sctl::Vector<Real>& center_in,
                                const sctl::Vector<Real>& radius_out,
                                const sctl::Vector<Real>& radius_in);





        // Accessors
        const sctl::Vector<Real>& centerCoordsOut() const { return _center_coords_out; }
        const sctl::Vector<Real>& centerCoordsIn() const { return _center_coords_in; }
        const sctl::Vector<Real>& radiusOut() const { return _radius_out; }
        const sctl::Vector<Real>& radiusIn() const { return _radius_in; }

};
#endif // __GLYMPH_WALL__