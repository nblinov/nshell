#pragma once

#include "nbody_system.hpp"

// Abstract integrator class which serves as a base for specific integrators like leapfrog
// Each specific implementation must have an update function that performs the timestep
class integrator
{
    public:
        nbody_system& gas; // Reference to the system of particles
        double t; // Current simulation time
        double dt;  // Current time-step

        // Initialize with the system of particles only
        integrator(nbody_system &gas_) : gas(gas_), t(0.), dt(0.) {} 
        // Initialize with system of particles, initial time and initial timestep
        integrator(nbody_system &gas_, const double t_, const double dt_) : gas(gas_), t(t_), dt(dt_) {} 
        
        // Perform time-step
        virtual void update()=0;
};

class leapfrog : public integrator
{
    public:
        // Initialize with the system of particles only
        leapfrog(nbody_system &gas_) : integrator(gas_) {}; 
        // Initialize with system of particles, initial time and initial timestep
        leapfrog(nbody_system &gas_, const double t_, const double dt_) : integrator(gas_, t_, dt_) {}; 

        void update()
        {
            // Leapfrog (technically this should be equivalent to verlet)
            // v(n+1/2) = v(n) + a(n)*dt/2
            // x(x+1) = x(n) + v(n+1/2)*dt
            // v(n+1) = v(n+1/2) + a(n+1)*dt/2
            // Only the integer step velocity and positions are stored (not half step)
            
            #pragma omp parallel for
            for (size_t n = 0; n < gas.size(); ++n)
            {
                gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
                gas[n].r += gas[n].vr*dt; // x(n) -> x(n+1)
                gas[n].a = gas.compute_acceleration(n); // a(n+1)
                gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 

                if (gas[n].r < 0.)
                  gas[n].r = 0.;
                gas[n].t += dt;
            }
            t += dt;
        }
};


