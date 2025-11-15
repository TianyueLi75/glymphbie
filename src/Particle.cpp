#include "Particle.hpp"

template <class Real>
Particle<Real>::Particle(const sctl::Vector<Real> X_, const Real T_, const sctl::Long Nt_) : T(T_), Nt(Nt_) {
    X.ReInit(Nt);
    X[0] = X_;
    dt = T / Nt;
    Nptcl = X_.Dim()/3;
    U.ReInit(Nt);
}

template <class Real>
Particle<Real>::~Particle() {} 

template <class Real> 
void Particle<Real>::SetParticleCoord(const sctl::Vector<Real> X_) {
    if (cur_tstep != 0) {
        std::cout << "Warning: attempting to change particle locations when time evolution is already underway. Aborted." << std::endl;
    } else {
        X[0] = X_;
        Nptcl = X_.Dim()/3;

    }
}

template <class Real>
void Particle<Real>::SetTime(const Real T_, const sctl::Long Nt_) {
    if (cur_tstep != 0) {
        std::cout << "Warning: attempting to change particle locations when time evolution is already underway. Aborted." << std::endl;
    } else {
        T = T_;
        Nt = Nt_;
        dt = T/Nt;
    }
}

template <class Real>
void Particle<Real>::EvolveInTime(sctl::Long time_solver, sctl::Long collision_mode, Annular<Real> channel, StokesBIO<Real> LayerPotenOp, sctl::Vector<Real> sigma) {
    SCTL_ASSERT(T > 0);
    SCTL_ASSERT(dt > 0);
    SCTL_ASSERT(cur_tstep == 0); 

    while (cur_tstep < Nt) {
        // Compute U at current particle position
        sctl::Vector<Real> Ucur(Nptcl*3, (sctl::Iterator<Real>) U[cur_tstep].begin(), false);
        LayerPotenOp.SetTargetCoord(X[cur_tstep]);
        LayerPotenOp.ComputePotential(&Ucur, sigma);

        cur_tstep += 1; // cur_tstep==0 is initial geom, don't change.
        
        if (time_solver == 0) {
            Euler_step(Ucur);
        } else {
            RK4_step(Ucur);
        }
        
        if (collision_mode==0) {
            // no collision
            continue;
        } else if (collision_mode == 1) {
            ResolveCollision1(channel);
        } else if (collision_mode == 2) {
            ResolveCollision2(channel);
        }
    }
}


template <class Real>
void Particle<Real>::Euler_step(sctl::Vector<Real> U) {
    sctl::Vector<Real> Xlast(Nptcl*3, (sctl::Iterator<Real>)X[cur_tstep-1].begin(),false);
    sctl::Vector<Real> Xcur(Nptcl*3, (sctl::Iterator<Real>)X[cur_tstep].begin(),false);
    Xcur = Xlast + dt * U;
}

template <class Real>
void Particle<Real>::RK4_step(sctl::Vector<Real> U) {
    // TODO
}

template <class Real>
void Particle<Real>::ResolveCollision1(Annular<Real> channel) {
    // No penetration: cuts off any particle outside channel to stop at boundary
    // check which points are outside domain
    sctl::Vector<Real> Xcur(Nptcl*3, (sctl::Iterator<Real>)X[cur_tstep].begin(), true);
    sctl::Vector<Real> Xupdated(Nptcl*3, (sctl::Iterator<Real>)X[cur_tstep].begin(), false);
    sctl::Vector<sctl::Long> in_channel = channel.InDomain(Xcur); // = 1 if in domain, = 0 otherwise
    sctl::Vector<Real> Xprev(Nptcl*3, (sctl::Iterator<Real>)X[cur_tstep-1].begin(), false); // values from previous step.
    Xupdated = Xcur * in_channel + Xprev * (1-in_channel);
}

template <class Real>
void Particle<Real>::ResolveCollision2(Annular<Real> channel)  {
    // opposite and damped repulsion force in normal direction.
    // TODO: need BIO, etc, or at least the U values.
}
