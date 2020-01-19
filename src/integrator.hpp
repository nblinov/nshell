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
            // Leapfrog 
            // v(n+1/2) = v(n) + a(n)*dt/2
            // x(x+1) = x(n) + v(n+1/2)*dt
            // v(n+1) = v(n+1/2) + a(n+1)*dt/2
            // Only the integer step velocity and positions are stored (not half step)
            
            // Keeps track which shells to update (and therefore which to skip)
            vector<size_t> shells_to_update;

            #pragma omp parallel for
            for (size_t n = 0; n < gas.size(); ++n)
            {
                // Time-step takes us to negative r 
                // Velocity changes sign when we go through the origin
                if (gas[n].r + (gas[n].vr+gas[n].a*dt/2)*dt <= 0.)
                {
                  gas[n].vr = -gas[n].vr;
                  gas[n].t += dt;
                  continue;
                }
                else
                {
                  // First leapfrog step using old accelerations
                  gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
                  gas[n].r += gas[n].vr*dt; // x(n) -> x(n+1)

                  shells_to_update.push_back(n);
                }
                /*
                // An alternative to the above: when vr -> -vr, continue time step as normal. Which one is right?
                // Time-step takes us to negative r 
                // Velocity changes sign when we go through the origin
                if (gas[n].r + (gas[n].vr+gas[n].a*dt/2)*dt <= 0.)
                {
                  gas[n].vr = -gas[n].vr;
                }
                // First leapfrog step using old accelerations
                gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
                gas[n].r += gas[n].vr*dt; // x(n) -> x(n+1)

                shells_to_update.push_back(n);
                */

            }

            // Now that all positions have been updated, compute new accelerations, velocities
            #pragma omp parallel for
            for (size_t i = 0; i < shells_to_update.size(); ++i)
            {
                  size_t n = shells_to_update[i];
                  gas[n].a = gas.compute_acceleration(n); // a(n+1)
                  gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 

                  gas[n].t += dt;
            }
            t += dt;
        }
};

// Same as above, but trying to vary the time step based on the smallest dynamical time among the shells
class adaptive_leapfrog : public integrator
{

    public:
        double c_dyn;  // Coefficient that determines the relation between dt and the smallest dynamical time in get_time_step
        double dt_min; // Minimum allowed time step
        double dt_max; // Maximum allowed time step
        double dt_prev; // previous time-step

        // Initialize with the system of particles only
        adaptive_leapfrog(nbody_system &gas_) : integrator(gas_) { initialize_timestep(); }; 
        // Initialize with system of particles, initial time and initial timestep
        adaptive_leapfrog(nbody_system &gas_, const double t_, const double dt_) : integrator(gas_, t_, dt_), c_dyn(0.01), dt_min(1e-5), dt_max(1e-1) { initialize_timestep(); }; 

        double get_time_step()
        {
          double smallest_t_dyn;
          smallest_t_dyn = gas[0].t_dyn;
          for(size_t i = 0; i < gas.size(); ++i)
          {
            if (gas[i].t_dyn < smallest_t_dyn)
              smallest_t_dyn = gas[i].t_dyn;
          }
          //cout << "smallest t_dyn = " << smallest_t_dyn << endl;
          return std::min(std::max(c_dyn*smallest_t_dyn, dt_min), dt_max);
        }

        void initialize_timestep()
        {
          dt = get_time_step(); 
        }

        void update()
        {
            // Leapfrog 
            // v(n+1/2) = v(n) + a(n)*dt/2
            // x(x+1) = x(n) + v(n+1/2)*dt
            // v(n+1) = v(n+1/2) + a(n+1)*dt/2
            // Only the integer step velocity and positions are stored (not half step)
            
            // Holder, Leimkuhler, Reich (2001) 
            //dt_prev = dt;
            dt = get_time_step(); // Naive
            //dt = 1./(2./get_time_step() - 1./dt_prev); // https://www.sciencedirect.com/science/article/abs/pii/S0168927401000897
            //dt = pow(get_time_step(),2.)/dt_prev; // https://arxiv.org/abs/1105.1082 
            //cout << "time step = " << get_time_step() << "   " << dt_prev << endl;

            // Keeps track which shells to update (and therefore which to skip)
            vector<size_t> shells_to_update;

            #pragma omp parallel for
            for (size_t n = 0; n < gas.size(); ++n)
            {
                // Time-step takes us to negative r 
                // Velocity changes sign when we go through the origin
                if (gas[n].r + (gas[n].vr+gas[n].a*dt/2)*dt <= 0.)
                {
                  gas[n].vr = -gas[n].vr;
                  gas[n].t += dt;
                  continue;
                }
                else
                {
                  // First leapfrog step using old accelerations
                  gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
                  gas[n].r += gas[n].vr*dt; // x(n) -> x(n+1)

                  shells_to_update.push_back(n);
                }
            }

            // Now that all positions have been updated, compute new accelerations, velocities
            #pragma omp parallel for
            for (size_t i = 0; i < shells_to_update.size(); ++i)
            {
                  size_t n = shells_to_update[i];
                  gas[n].a = gas.compute_acceleration(n); // a(n+1)
                  gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 

                  gas[n].t += dt;
            }
            t += dt;

        }
};

