#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <csbq.hpp>

// TODO: time solver mode / collision mode use enumerate instead of numbers?
// TODO: add Comm.

template <class Real> class Particle {

    public:

        Particle(const sctl::Vector<Real> X_, const Real T_, const sctl::Long Nt_);
        ~Particle();

        void SetParticleCoord(const sctl::Vector<Real> X_);
        void SetTime(const Real T_, const sctl::Long Nt_);
        void EvolveInTime(sctl::Long time_solver, sctl::Long collision_mode, Annular<Real> channel, StokesBIO<Real> LayerPotenOp, sctl::Vector<Real> sigma);

    private:
        sctl::Vector<sctl::Vector<Real>> X; // vector of vectors of particle positions at each time step
        sctl::Vector<sctl::Vector<Real>> U;
        sctl::Long Nptcl;
        Real dt = -1.; // Time step of simulation; will always be initialized by constructor.
        Real T = -1.; // Total time of simulation
        sctl::Long Nt;
        sctl::Long cur_tstep = 0;

        void Euler_step(sctl::Vector<Real> U);
        void RK4_step(sctl::Vector<Real> U);
        void ResolveCollision1(Annular<Real> channel);
        void ResolveCollision2(Annular<Real> channel);

}
#endif