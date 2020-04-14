#pragma once

#include "nbody_system.hpp"

// This class contains all the particles of the simulation and defines 
// various functions to compute forces and other properties for the EMD scenario
// EMD particle lifetime is specified by tau_emd; the fraction of stable matter is stable_frac
// This class inherits from nbody_system so we only have to redefine the computation 
// of the background energy density in radiation and the calculation of the shell masses 
// which now depend on time.
class nbody_system_emd : public nbody_system
{
    public:
    const double tau_emd; // Lifetime of the particle respeonsible for EMD
    const double stable_frac; // Fraction of DM that is stable (that does not decay when EMD ends)

    nbody_system_emd(size_t nshells_, const double tau_, const double sf_) : nbody_system(nshells_), tau_emd(tau_), stable_frac(sf_){} 
    nbody_system_emd(size_t nshells_, const double tau_, const double sf_, const double fs_) : nbody_system(nshells_,fs_), tau_emd(tau_), stable_frac(sf_){} 

    // Background density at time t in matter domination
    double get_density_bg(double t)
    {
      //cout << "Using EMD bg density!" << endl;
      return 3./(2.*PI*t*(15.*tau_emd + 16.*t));
    }

    // Mass interior to radius r
    double get_mass_interior_to_r(double r)
    {
      double t = gas[0].t;  // Assume all particles have the same time coordinate
      double mass_interior = gas[0].mass; // Add a point mass at r = 0
      double rm, rp; // Shell boundaries
      // Finite shell thickness
      if (delta_r > 0.){
        update_shell_order();
        for (size_t m = 0; m < gas.size(); ++m)
        {
          rm = (m > 0) ? gas[shell_order[m-1]].r : 0.; // If inner-most shell, set left edge to r=0
          rp = (m < gas.size()-1) ? gas[shell_order[m+1]].r : gas[shell_order[m]].r; // If outer-most shell set right edge to shell position

          mass_interior += gas[m].mass*shell_mass_enclosed(r, rm, rp); // Shell thickness is r(i+1) - r(i-1)
        }
      }
      // Infinitely thin shells
      else {
        //return mass_interior + (n)*gas[0].mass;
        for (size_t m = 0; m < gas.size(); ++m)
          mass_interior += gas[m].mass*heaviside(r - gas[m].r);
      }

      //cout << "Using the EMD mass definition!" << endl;
      return mass_interior*(stable_frac + (1.-stable_frac)*exp(-t/tau_emd));
    }

};
