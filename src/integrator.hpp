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

            //#pragma omp parallel for default(none) shared(gas, shells_to_update, dt)
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
            //#pragma omp parallel for default(none) shared(gas, shells_to_update, dt)
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

// The methods below aren't working yet.
// Reversible adaptive scheme from https://www.sciencedirect.com/science/article/abs/pii/S0168927401000897 
// or using an alternative time step recursion: https://link.springer.com/content/pdf/10.1023/A:1022313123291.pdf
class reversible_adaptive_leapfrog : public integrator
{

    public:
        double ds;  // Coefficient that determines the relation between dt and the auxiliary parameter rho (see below) 
        double dt_min; // Minimum allowed time step
        double dt_max; // Maximum allowed time step
        double dt_prev; // previous time-step

        double rho; // Auxiliary variable from which we derive the current step
        double rho_prev;



        // Initialize with system of particles, initial time and initial timestep
        reversible_adaptive_leapfrog(nbody_system &gas_, const double t_, const double dt_) : integrator(gas_, t_, dt_), ds(0.01), dt_min(1e-5), dt_max(1e-1) { initialize_auxiliary(); }; 

        double scaling_function()
        {
          double smallest_t_dyn;
          smallest_t_dyn = gas[0].t_dyn;
          for(size_t i = 0; i < gas.size(); ++i)
          {
            if (gas[i].t_dyn < smallest_t_dyn)
              smallest_t_dyn = gas[i].t_dyn;
          }
          return 1./smallest_t_dyn;
          /*
          double sum;
          sum = 0.;
          for(size_t i = 0; i < gas.size(); ++i)
          {
            sum += gas[i].a*gas[i].a;
          }
          return sqrt(sum);
          */
        }

        void initialize_auxiliary()
        {
          //rho = scaling_function(); 
          rho = 1./scaling_function(); 
        }

        void update()
        {
            // Modified Leapfrog 
            // v(n+1/2) = v(n) + 0.5*a(n)*ds/rho(n)
            // x(x+1/2) = x(n) + 0.5*v(n+1/2)*ds/rho(n)
            // rho(n+1) = 2U(x(n+1/2)) - rho(n)
            // v(n+1) = v(n+1/2) + 0.5*a(n+1)*ds/rho(n+1)
            // x(n+1) = x(n+1/2) + 0.5*v(n+1/2)*ds/rho(n+1)
            
            //dt = ds/rho; 
            dt = ds*rho; 

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
                }
                
                // First leapfrog step using old accelerations (slightly different from normal leapfrog!)
                gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
                gas[n].r += gas[n].vr*dt/2; // x(n) -> x(n+1/2)

                shells_to_update.push_back(n);
              
            }

            //cout << rho << endl;
            rho_prev = rho;
            //rho = 2.*scaling_function() - rho_prev; // rho(n) -> rho(n+1)
            // dt = ds/rho;
            rho = 1./(2./scaling_function() - 1./rho_prev);
            dt = ds*rho;

            // Now that all positions have been updated, compute new accelerations, velocities
            #pragma omp parallel for
            for (size_t i = 0; i < shells_to_update.size(); ++i)
            {
                  size_t n = shells_to_update[i];
                  gas[n].r += gas[n].vr*dt/2; // x(n+1/2) -> x(n+1) 
            }

            #pragma omp parallel for
            for (size_t i = 0; i < shells_to_update.size(); ++i)
            {
                  size_t n = shells_to_update[i];

                  gas[n].a = gas.compute_acceleration(n); // a(n+1)
                  gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 
                  //gas[n].t += 0.5*ds*(1./rho + 1./rho_prev);
                  gas[n].t += 0.5*ds*(rho + rho_prev);
            }

            //t += 0.5*ds*(1./rho + 1./rho_prev);
            t += 0.5*ds*(rho + rho_prev);
            cout << "t = " << rho << " " << rho_prev << endl;

        }
};
