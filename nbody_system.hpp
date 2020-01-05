#pragma once

#include <cmath>
#include <cstddef>
using std::size_t;
#include <vector>
using std::vector;
#include "shell.hpp"
#include "utility.hpp"

class nbody_system
{
    public:
    size_t nshells;
    vector<shell> gas;
    const double force_softening; // units of distance
    const double force_softening_sq;

    nbody_system(size_t nshells_) : nshells(nshells_), gas(nshells_, shell(0.,0.,0.,0.)), force_softening(0.),  force_softening_sq(0.){} 
    nbody_system(size_t nshells_, const double fs_) : nshells(nshells_), gas(nshells_, shell(0.,0.,0.,0.)), force_softening(fs_),  force_softening_sq(fs_*fs_){} 

    shell& operator[](size_t i) { return gas[i]; }
    const shell& operator[](size_t i) const { return gas[i]; }

    size_t size()
    {
        return gas.size();
    }

    // Recomputes the acceleration and updates the dynamical time for the shell n
    double compute_acceleration(size_t n)
    {
        double a = 0.;
        double t, r;
        double mass_interior = 0.;

        t = gas[n].t;
        r = gas[n].r;

        mass_interior = get_mass_interior(n);
        gas[n].t_dyn = sqrt(5.*pow(r,3.)/(2.*mass_interior));

        a += + pow(gas[n].l,2.)/pow(r,3.) - mass_interior*r/pow(r*r+force_softening_sq,3./2.);
        return a;
    }
    // Background mass interior to r at time t in matter domination
    double get_mass_interior_bg(double t, double r)
    {
      return (2./9.)*pow(r,3)/pow(t,2);
    }

    // Mass interior to shell n
    double get_mass_interior(size_t n)
    {
      //double mass_interior = 0.;
      double mass_interior = gas[0].mass; // Add a point mass at r = 0
      //return mass_interior + (n)*gas[0].mass;
      for (size_t m = 0; m < gas.size(); ++m)
            mass_interior += gas[m].mass*heaviside(gas[n].r - gas[m].r);

      return mass_interior;
      /*
      */
    }

    // Mass interior to radius r
    double get_mass_interior_to_r(double r)
    {
      double mass_interior = 0.;
      for (size_t m = 0; m < gas.size(); ++m)
            mass_interior += gas[m].mass*heaviside(r - gas[m].r);

      return mass_interior;
    }

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

    inline double get_shell_energy(size_t n)
    {
      return 0.5*gas[n].vr*gas[n].vr + gas[n].l*gas[n].l/(gas[n].r*gas[n].r) - get_mass_interior(n)/gas[n].r;
    }

    double get_total_energy()
    {
      double total_energy;
      total_energy = 0.;
      for (size_t i = 0; i < gas.size(); i++)
        total_energy += get_shell_energy(i);

      return total_energy;
    }
};
