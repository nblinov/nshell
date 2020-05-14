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

#include <sstream>

#include "src/shell.hpp"
#include "src/nbody_system_emd.hpp"
#include "src/integrator.hpp"

int sim_steps = 0;                  // Integral simulation time. Real time given by dt*sim_time
const int nshells = 10000;  // Number of shell sot simulate


double epsilon = 0.5; // The initial perturbation has profile delta ~ (M/M_0)^-epsilon

// Parameter that specifies angular momentum of the shells following https://arxiv.org/pdf/astro-ph/0008217.pdf
double alpha = 0.01;  // if alpha < 1, then IC is such that shell is bound

double delta_avg_init = 1e-1; // Initial average overdensity of the whole shell config (out to radius rmax)
const double Mtot = 1.; // Arbitrary normalization of shell masses
const double rmax = 1.; // Arbitrary normalization of distances

// Infer start time from initial overdensity
const double dt = 0.00005; // Initial stepsize (may be modified in the integrator)
const double tau_emd_over_ti = 100.;

// These need to be updated if delta_avg_init is changed 
double ti = sqrt((2./9.) * (1.+delta_avg_init) * (1. + 1./sqrt(tau_emd_over_ti)) * exp(-1./tau_emd_over_ti)); 
double tau_emd = tau_emd_over_ti*ti; // lifetime of the field responsible for EMD 
double sim_time_max = 50.*tau_emd;  // Length of the simulation

const double stable_frac = 0.0001; // fraction of DM that does NOT decay at the end of EMD

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
    double mass_interior;
    double delta;

	for(size_t i = 0; i < nshells; i++){

      mass_interior = (i+1)*m*(stable_frac + (1.-stable_frac)*exp(-1/tau_emd_over_ti)); 
      delta = delta_avg_init/pow((i+1.)/nshells,epsilon);

      ri = pow(mass_interior*(1. + delta_avg_init)/(1. + delta),1./3.); 
      vri = (2./3.)*(ri/ti)*(1 - 1./3.);
       

      // Initialize angular momentum
      l = sqrt(alpha*2.*ri*mass_interior*delta);

      gas[i] = shell(m, ri, vri, l, ti);

	}

    // Which shells to print out to standard output via output_shell_evolution below
    r_out.push_back(1);
    r_out.push_back(nshells/2);
    r_out.push_back(nshells-1);

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

    if (argc < 3){
      cout << "Need three arguments:  output directory, initial profile slope epsilon and angular momentum parameter alpha!" << endl;
      exit(0);
    }

    std::string out_dir = argv[1];  
    epsilon = atof(argv[2]);
    alpha = atof(argv[3]);

    // Set the initial overdensity, and re-compute times
    if (argc > 3){
      delta_avg_init = atof(argv[4]);
      ti = sqrt((2./9.) * (1.+delta_avg_init) * (1. + 1./sqrt(tau_emd_over_ti)) * exp(-1./tau_emd_over_ti)); 
      tau_emd = tau_emd_over_ti*ti; // lifetime of the field responsible for EMD 
      sim_time_max = 50.*tau_emd;  // Length of the simulation
    }

    std::string out_param = out_dir + "/log.param";

    // Write down the parameters of the run
    ofstream param_file (out_param.c_str());
    param_file << "Initial profile slope epsilon = " << epsilon << endl;
    param_file << "Angular momentum parameter alpha = " << alpha << endl;
    param_file << "Number of shells = " << nshells << endl;
    param_file << "Initial average overdensity = " << delta_avg_init << endl;
    param_file << "Stable DM fraction = " << stable_frac << endl;
    param_file << "Minimum time step (if adaptive) = " << dt << endl;

    // Initialize the n-shell system and specify initial conditions
    nbody_system_emd gas(nshells, tau_emd, stable_frac);
    initialize_gas(gas);

    // Some sanity checks
    param_file << "Sanity checks..." << endl;
    param_file << "Simulating " << gas.size() << " shells..." << endl;
    param_file << "Max radius = " << gas[nshells-1].r << "; M(nshells) = " << gas.get_mass_interior(nshells-1) << endl;
    param_file << "Average overdensity = " << gas.get_delta_interior(ti, nshells-1) << "; should be = " << delta_avg_init << endl;

    param_file.close();


    // Choose the integrator
    //leapfrog stepper(gas, ti, dt); 
    adaptive_leapfrog stepper(gas, ti, dt); 

    // Run the simulation and print output
    while(stepper.t < sim_time_max)
    {
        stepper.update();

        // Output a few individual shell properties
        if ((sim_steps%200 == 0)) output_shell_evolution(gas);

        // Output binned density profile of the whole system
        //if ((sim_steps%100000 == 0)) 
        if ((sim_steps%1000 == 0)) 
        {
          std::string out_name;
          std::ostringstream ss;
          //ss << floor(stepper.t); 
          ss << round(100.*gas[0].t/tau_emd)/100.;
          out_name = out_dir+"/radial_profile_"+ss.str()+".dat"; 
          gas.output_radial_profile(out_name.c_str());
        }

        sim_steps++;

    }
}
