#pragma once

#include <cmath>
#include <cstddef>
using std::size_t;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;

#include <algorithm>
using std::sort;

#include "shell.hpp"
#include "utility.hpp"


// This class contains all the particles of the simulation and defines 
// various functions to compute forces and other properties
class nbody_system
{
    // Functor used in radius-ordering shells
    struct shell_comparator{
      nbody_system& sys;
      shell_comparator(nbody_system& sys_) : sys(sys_) {}
      bool operator() (size_t i, size_t j) { return (sys[i].r < sys[j].r);}
    };

    public:
    size_t nshells; // Number of shells in the system
    vector<shell> gas;  // Container for the particles;
    vector<size_t> shell_order; // Vector of indices for particles in gas, ordered according to shell radius
    const double force_softening; // Force softening parameter (units of distance)
    const double force_softening_sq; // Square of the above

    double delta_r; // Shell thickness (if > 0, take shells to be equal density between neighbouring shells; actual size of delta_r ignored)

    nbody_system(size_t nshells_) : nshells(nshells_), gas(nshells_, shell(0.,0.,0.,0.)), 
                                    force_softening(0.),  force_softening_sq(0.), shell_order(nshells_,0), delta_r(0.){} 
    nbody_system(size_t nshells_, const double fs_) : nshells(nshells_), gas(nshells_, shell(0.,0.,0.,0.)), 
                                                      force_softening(fs_),  force_softening_sq(fs_*fs_), shell_order(nshells_,0), delta_r(0.){} 
    // Access individual particles
    shell& operator[](size_t i) { return gas[i]; }
    const shell& operator[](size_t i) const { return gas[i]; }

    size_t size()
    {
        return gas.size();
    }

    void set_shell_thickness(double dr)
    {
      delta_r = dr;
    }

    void update_shell_order()
    {
      sort(shell_order.begin(), shell_order.end(), shell_comparator(*this));
    }

    // Computes the acceleration and updates the dynamical time for the shell n
    virtual double compute_acceleration(size_t n)
    {
        double a = 0.;
        double t, r;
        double mass_interior = 0.;

        t = gas[n].t;
        r = gas[n].r;

        mass_interior = get_mass_interior(n);
        gas[n].t_dyn = sqrt(5.*pow(r,3.)/(2.*mass_interior));
        //cout << "t_dyn = " << gas[n].t_dyn << "   "  << r << "   " << mass_interior << endl;

        a += + pow(gas[n].l,2.)/pow(r,3.) - mass_interior*r/pow(r*r+force_softening_sq,3./2.);//- r;
        return a;
    }

    // Background density at time t in matter domination
    virtual double get_density_bg(double t)
    {
      return 1./(6.*PI*t*t);
    }

    // Background mass interior to r
    virtual double get_mass_interior_bg(double t, double r)
    {
      return (4.*PI/3.)*pow(r,3.)*get_density_bg(t); //(2./9.)*pow(r,3)/pow(t,2);
    }

    // Fraction of shell (shell boundaries are between rm and rp, rp > rm) enclosed at radius r
    double shell_mass_enclosed(double r, double rm, double rp)
    {
      double norm = pow(rp,3.) - pow(rm,3.);
      return std::max(std::min((pow(r,3.) - pow(rm,3.))/norm,1.),0.);
    }

    // Mass interior to radius r
    virtual double get_mass_interior_to_r(double r)
    {
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
          //mass_interior += gas[m].mass*shell_mass_enclosed(r, gas[m].r-delta_r/2., gas[m].r+delta_r/2.); // Shell thickness is fixed to be delta_r
          //cout << shell_mass_enclosed(r, rm, rp) << endl;
        }
      }
      // Infinitely thin shells
      else {
        //return mass_interior + (n)*gas[0].mass;
        for (size_t m = 0; m < gas.size(); ++m)
          mass_interior += gas[m].mass*heaviside(r - gas[m].r);
      }

      return mass_interior;
    }

    // Mass interior to shell n
    double get_mass_interior(size_t n)
    {
      return get_mass_interior_to_r(gas[n].r);
    }

    // Compute density at r via a discrete derivative
    double get_density(double r, double dr)
    {
      double r1 = r - dr/2;
      double r2 = r + dr/2;

      double mass_interior1 = get_mass_interior_to_r(r1);
      double mass_interior2 = get_mass_interior_to_r(r2);

      return (mass_interior2-mass_interior1)/dr/(4.*PI*r*r);
    }

    // Overdensity interior to shell n at time t
    inline double get_delta_interior(double t, size_t n)
    {
      return get_mass_interior(n)/get_mass_interior_bg(t, gas[n].r) - 1.;
    }

    // Overdensity interior to radius r at time t
    inline double get_delta_interior_to_r(double t, double r)
    {
      return get_mass_interior_to_r(r)/get_mass_interior_bg(t, r) - 1.;
    }

    inline double get_shell_energy(size_t n)
    {
      return 0.5*gas[n].vr*gas[n].vr + gas[n].l*gas[n].l/(gas[n].r*gas[n].r) - get_mass_interior(n)/sqrt(gas[n].r*gas[n].r + force_softening_sq); //+ 0.5*gas[n].r*gas[n].r;//;
    }

    double get_total_energy()
    {
      double total_energy;
      total_energy = 0.;
      for (size_t i = 0; i < gas.size(); i++)
        total_energy += get_shell_energy(i);

      return total_energy;
    }

    double get_r_max()
    {
      double largest_r;
      largest_r = gas[0].r;
      for(size_t i = 0; i < nshells; ++i)
      {
        if (gas[i].r > largest_r)
          largest_r = gas[i].r;
      }
      return largest_r; 
    }

    virtual void output_radial_profile(const char * file_name)
    {
        double r, dr;
        int nbins;
        nbins = 20;
        dr = 1.1*get_r_max()/nbins;
        ofstream radial_profile (file_name);
        for (size_t i = 0; i < nbins; i++)
        {
          r = (i+1)*dr;
          radial_profile << r << "    " << get_density(r, dr)/get_density_bg(gas[0].t) << "    " << get_mass_interior_to_r(r) << "    " << get_delta_interior_to_r(gas[0].t,r) << endl;  
        }
        radial_profile.close();
    }
};
