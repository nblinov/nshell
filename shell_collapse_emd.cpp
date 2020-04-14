/*
 * In this example we evaluate the collapse of a constant density sphere to the turn-around point 
 * during EMD followed by RD 
 * The initial radial velocity of the shells is given by vr = H*r
 * This turn-around point (radius and time) agree with those estimated based on analytics for these same initial conditions
 */
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cstdlib>	
using std::exit;

#include <vector>
using std::vector;

#include "src/shell.hpp"
#include "src/nbody_system_emd.hpp"
#include "src/integrator.hpp"

int sim_steps = 0;                  // Integral simulation time. Real time given by dt*sim_time
const double sim_time_max = 2000.;  // Length of the simulation
const int nshells = 100;  // Number of shell sot simulate


const double one_plus_delta_init = 1. + 1e-2; // Initial overdensity of the whole shell config
const double Mtot = 1.; // Arbitrary normalization of shell masses
const double rmax = 1.; // Arbitrary normalization of distances
const double Mtot_interior = Mtot*(nshells-1.)/nshells;

// Infer start time from initial overdensity
const double ti = sqrt((2./9.) * one_plus_delta_init * pow(rmax,3.)/Mtot);
const double dt = 0.0005; // Initial stepsize (may be modified in the integrator)
const double tau_emd = 100.; // lifetime of the field responsible for EMD 
const double stable_frac = 0.01; // fraction of DM that does NOT decay at the end of EMD

vector<double> r_avg; // Radii for which to print out 1+delta(r)
vector<size_t> r_out; // Indices of radii for which to print out 1+delta(r)

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
    double alpha = 0.01;  // if alpha < 1, then IC is such that shell is bound
    double mass_interior;

	for(size_t i = 0; i < nshells; i++){
      //ri = (i+1.)/(nshells+1.)*rmax; // Equally spaced initial conditions: impossible to make all shells start out in the linear regime

      // These initial conditions assume constant radial profile for the overdensity
      // i.e. delta is the same at every radius
      ri = pow((9./2.)*ti*ti*(i+1)*m/one_plus_delta_init,1./3.); 
      vri = (2./3.)*(ri/ti)*(1 - 0*(one_plus_delta_init-1.)/3.);
       
      // will update the angular momentum in the nex loop
      l = ri;

      gas[i] = shell(m, ri, vri, l, ti);

      // Initialize angular momentum
      mass_interior = gas.get_mass_interior(i);
      l = sqrt(alpha*2.*gas[i].r*mass_interior*(one_plus_delta_init-1.));
      gas[i].l = l;

	}

    //cout << "At ti, mass interior to rmax, bg = " << ti << "\t" << get_mass_interior_to_r(rmax)/get_mass_interior_bg(ti, rmax) << endl;
	cout << "# Simulating " << gas.size() << " shells..." << endl;

    r_out.push_back(1);
    r_out.push_back(nshells/2);
    r_out.push_back(nshells-1);

    cout << "# Expected turn-around times and radii for shells are " << endl;
    for (size_t m = 0; m < r_out.size(); m++){
      double delta_i = gas.get_delta_interior(ti, r_out[m]);
      // These turn-around estimates are based on radial collapse, with IC corresponding to v_r = Hr
      cout <<  "# ri = " << r_out[m] << " delta = " <<  delta_i << " r_ta = " << gas[r_out[m]].r/(delta_i) << " t_ta = " << ti*3*PI/(4.*pow(delta_i, 3./2.)) << endl;
      cout <<  "# vr = " << gas[r_out[m]].vr << " vtheta = " << gas[r_out[m]].l/gas[r_out[m]].r << endl;
    }
    cout << "# Max radius = " << gas[nshells-1].r << " contains M = " << gas.get_mass_interior(nshells-1) << endl;
}

void output_shell_evolution(nbody_system &gas)
{
  size_t i;
  cout << std::scientific; 
  double t = gas[0].t;

  cout << t << "    ";

  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    cout << gas[i].r << "    " << gas[i].vr << "    " << gas.get_delta_interior(t, i) << "    " << gas.get_mass_interior(i) << "    ";

  }
  //cout << get_total_energy() << "    ";
  double t_dyn;
  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    t_dyn = gas[i].t_dyn;
    cout << t_dyn/dt << "    ";
  }

  cout << endl;
}

int main(int argc, char **argv)
{
    // Initialize the n-shell system and specify initial conditions
    nbody_system_emd gas(nshells, tau_emd, stable_frac);
	initialize_gas(gas);

    // Choose the integrator
    leapfrog stepper(gas, ti, dt); 

    // Run the simulation and print output
    while(stepper.t < sim_time_max)
    {
        stepper.update();
        if ((sim_steps%100 == 0)) output_shell_evolution(gas);
        sim_steps++;
    }
}
