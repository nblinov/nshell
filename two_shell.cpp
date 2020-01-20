/*
 * In this example we look at the evolution of two shells 
 */
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cstdlib>	
using std::exit;

#include <vector>
using std::vector;

#include <string>
#include <sstream>


#include "src/shell.hpp"
#include "src/nbody_system.hpp"
#include "src/integrator.hpp"

// Simulation parameters
int sim_steps = 0;                  // Integral simulation time. Real time given by dt*sim_time
const double sim_time_max = 5000.;  // Length of the simulation
const int nshells = 2;  // Number of shells to simulate

const double epsilon = 0.5; // The initial perturbation has profile delta ~ (M/M_0)^-epsilon
const double delta_at_r0 = 1e-1; // Initial overdensity at the origin
const double Mtot = 1.; // Arbitrary normalization of shell masses
const double rmax = 1.; // Arbitrary normalization of distances
const double Mtot_interior = Mtot*(nshells-1.)/nshells;

// Infer start time from initial overdensity
const double ti = sqrt((2./9.) * (1. + delta_at_r0) * pow(rmax,3.)/Mtot);
const double dt = 0.001; // Initial stepsize (may be modified in the integrator)

// Output-related variables
vector<double> r_avg; // Radii for which to print out 1+delta(r)
vector<size_t> r_out; // Indices of radii for which to print out detailed output
vector<double> r_ta;  // Turn-around radii for the above
vector<double> t_ta;  // Turn-around times for the above


// This function sets up the initial conditions
void initialize_gas(nbody_system &gas)
{
  double m;
  double ri;
  double vri;
  double l;

// Shell initialization:
// shell(mass, double position, double velocity, double angular_momentum, double initial_time);

    // We normalize shell masses such that the total mass is Mtot
    m = Mtot/(nshells);

    // Initializing following https://arxiv.org/pdf/astro-ph/0008217.pdf
    double alpha = 0.0;  // Parametrization of shell angular momentum; if alpha < 1, then IC is such that shell is bound
    double mass_interior;
    double delta;

	for(size_t i = 0; i < nshells; i++){
      mass_interior = (i+1)*m;
      delta = delta_at_r0/pow(i+1.,epsilon);
      
      // initial shell radius
      ri = pow((9./2.)*ti*ti*mass_interior/(1.+delta),1./3.); 
      // initial shell radial velocity
      vri = (2./3.)*(ri/ti);
      // conserved shell angular momentum
      l = sqrt(alpha*2.*ri*mass_interior*delta);

      gas[i] = shell(m, ri, vri, l, ti);
      gas[i].a = gas.compute_acceleration(i);
      gas.shell_order[i] = i;

      r_out.push_back(i);
	}

    // Thick shells
    //gas.set_shell_thickness(1.0); //dr > 0 => thick shells
    //gas.set_shell_thickness(gas[nshells-1].r/nshells);

    //cout << "At ti, mass interior to rmax, bg = " << ti << "\t" << get_mass_interior_to_r(rmax)/get_mass_interior_bg(ti, rmax) << endl;
	cout << "# Simulating " << gas.size() << " shells..." << endl;

    cout << "# Expected turn-around times and radii for shells are " << endl;
    for (size_t m = 0; m < r_out.size(); m++){
      double delta_i = gas.get_delta_interior(ti, r_out[m]);

      cout <<  "# ri = " << r_out[m] << "; vr = " << gas[r_out[m]].vr << "; vtheta = " << gas[r_out[m]].l/gas[r_out[m]].r << "; a = " <<  gas[r_out[m]].a << endl;
    }
    cout << "# Max radius = " << gas[nshells-1].r << " contains M = " << gas.get_mass_interior(nshells-1) << endl;
    //exit(0);
}

// Prints times, radii, velocities and enclosed masses and overdensities for a few shells
void output_shell_evolution(nbody_system& gas, integrator& stepper)
{
  size_t i;
  cout << std::scientific; 
  double t;
  double current_dt;

  // Current timestep
  current_dt = stepper.dt;

  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    t = gas[i].t;
    cout << t << "    " << gas[i].r << "    " << gas[i].vr << "    " << gas[i].a << "    ";
    //cout << gas[i].r << "    " << gas[i].vr << "    " << gas.get_delta_interior(t, i) << "    " << gas.get_mass_interior(i) << "    ";
  }
  double t_dyn;
  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    t_dyn = gas[i].t_dyn;
    cout << t_dyn/current_dt << "    ";
  }
  cout << gas.get_total_energy() << "    ";
  cout << current_dt;

  cout << endl;
}

int main(int argc, char **argv)
{
    const double force_softening = .1;

    // Initialize the n-shell system and specify initial conditions
    nbody_system gas(nshells, force_softening);
	initialize_gas(gas);

    // Choose the integrator
    //leapfrog stepper(gas, ti, dt); 
    adaptive_leapfrog stepper(gas, ti, dt); 
    // This one is not working yet
    //reversible_adaptive_leapfrog stepper(gas, ti, dt); 

    // Run the simulation and print output
    while(stepper.t < sim_time_max)
    {
        stepper.update();

        // Output evolution for a few shells every few steps
        if ((sim_steps%100 == 0)) output_shell_evolution(gas, stepper);

        sim_steps++;
    }
}
