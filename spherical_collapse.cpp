#include <cstddef>
using std::size_t;

#include <cassert>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ofstream;

#include <iomanip>
using std::setw;

#include <cmath>

#include <cstdlib>	
using std::exit;

#include <vector>
using std::vector;

#include "shell.hpp"

#include "mtrand.hpp"
mtrand Rand;

vector<shell> gas; // gas = vector of particles

const int nshells = 100;  // Number of shell sot simulate
const double force_softening_alpha = 0.05; // units of distance
const double force_softening_alpha2 = force_softening_alpha*force_softening_alpha;


double dt = 0.0005;
double dt2 = dt*dt;
const double c_dyn = 0.01;  // Coefficient that determines the relation between dt and the smallest dynamical time in get_time_step
const double dt_min = 1e-5; // Minimum allowed time step
const double dt_max = 1e-1; // Maximum allowed time step

const double one_plus_delta_init = 1. + 1e-2; // Initial overdensity of the whole shell config

const double Mtot = 1.; // Arbitrary normalization of shell masses
const double rmax = 1.; // Arbitrary normalization of distances
const double Mtot_interior = Mtot*(nshells-1.)/nshells;

// Infer start time from initial overdensity
const double ti = sqrt((2./9.) * one_plus_delta_init * pow(rmax,3.)/Mtot);


const int delay = 0.1; // milliseconds


vector<double> r_avg; // Radii for which to print out 1+delta(r)
vector<size_t> r_out; // Indices of radii for which to print out 1+delta(r)

int sim_time = 0;				 // Integral simulation time. Real time given by dt*sim_time
const double sim_time_max = 5000; 		 // Length of the simulation

#ifdef USE_OPEN_GL
void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	for (size_t n = 0; n < gas.size(); ++n)
		gas[n].draw();
	glutSwapBuffers();
}
#endif // USE_OPEN_GL	

inline double heaviside(double x){
  return 1. ? x > 0. : 0.;
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

// Overdensity interior to shell n at time t
inline double get_delta_interior(double t, size_t n)
{
  return get_mass_interior(n)/get_mass_interior_bg(t, gas[n].r) - 1.;
}

void initialize_gas(void)
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
    double alpha = 0.5;  // if alpha < 1, then IC is such that shell is bound
    double mass_interior;

    alpha = 0.1;

	for(size_t i = 0; i < nshells; i++){
      //ri = (i+1.)/(nshells+1.)*rmax; // Equally spaced initial conditions: impossible to make all shells start out in the linear regime

      // These initial conditions assume constant radial profile for the overdensity
      // i.e. delta is the same at every radius
      ri = pow((9./2.)*ti*ti*(i+1)*m/one_plus_delta_init,1./3.); 
      vri = (2./3.)*(ri/ti)*(1 - (one_plus_delta_init-1.)/3.);
       
      // will update the angular momentum in the nex loop
      l = ri;

      //cout << i << "\t" << ri << "\t" << vri << endl;

      gas.push_back(shell(m, ri, vri, l, ti));

      // Initialize angular momentum
      mass_interior = get_mass_interior(i);
      l = sqrt(alpha*2.*gas[i].r*mass_interior*(one_plus_delta_init-1.));
      gas[i].l = l;

	}

    //cout << "At ti, mass interior to rmax, bg = " << ti << "\t" << get_mass_interior_to_r(rmax)/get_mass_interior_bg(ti, rmax) << endl;
	cout << "# Simulating " << gas.size() << " shells..." << endl;

    // r_avg.push_back(0.1*rmax);
    // r_avg.push_back(0.5*rmax);
    // r_avg.push_back(rmax);
    r_out.push_back(1);
    r_out.push_back(nshells/2);
    r_out.push_back(nshells-1);

    cout << "# Expected turn-around times and radii for shells are " << endl;
    for (size_t m = 0; m < r_out.size(); m++){
      double delta_i = get_delta_interior(ti, r_out[m]);
      // These turn-around estimates are based on radial collapse, with IC corresponding to v_r = Hr
      //cout <<  "# ri = " << r_out[m] << " delta = " <<  delta_i << " r_ta = " << gas[r_out[m]].r/(delta_i) << " t_ta = " << ti*3*3.14159/(4.*pow(delta_i, 3./2.)) << endl;
      // These turn-around estimates are based on radial collapse, with IC corresponding to v_r = Hr(1-delta/3) (from linear PT)
      cout <<  "# ri = " << r_out[m] << " delta = " <<  delta_i << " r_ta = " << 3.*gas[r_out[m]].r/(5.*delta_i) << " t_ta = " << ti*1.0951/pow(delta_i, 3./2.) << endl;
      cout <<  "# vr = " << gas[r_out[m]].vr << " vtheta = " << gas[r_out[m]].l/gas[r_out[m]].r << endl;
    }
    cout << "# Max radius = " << gas[nshells-1].r << " contains M = " << get_mass_interior(nshells-1) << endl;




}


// Recomputes the acceleration and updates the dynamical time for the shell n
double recompute_acceleration(size_t n)
{
	double a = 0.;
    double t, r;
    double mass_interior = 0.;

    t = gas[n].t;
    r = gas[n].r;

    mass_interior = get_mass_interior(n);
    gas[n].t_dyn = sqrt(5.*pow(r,3.)/(2.*mass_interior));

    a += + pow(gas[n].l,2.)/pow(r,3.) - mass_interior*r/pow(r*r+force_softening_alpha2,3./2.);
	return a;
}

double get_time_step()
{
  double smallest_t_dyn;
  smallest_t_dyn = gas[0].t_dyn;
  for(size_t i = 0; i < nshells; ++i)
  {
    if (gas[i].t_dyn < smallest_t_dyn)
      smallest_t_dyn = gas[i].t_dyn;
  }
  return std::min(std::max(c_dyn*smallest_t_dyn, dt_min), dt_max);
}
/*
*/

void verlet_update()
{
	// self-starting Verlet algorithm ~ Leap-frog algo
	// x(n+1) = x(n) + v(n)dt + 0.5*a(n)dt^2
	// v(n+1) = v(n) + 0.5*(a(n+1)+a(n))*dt
    #pragma omp parallel for
	for (size_t n = 0; n < gas.size(); ++n)
    {
		gas[n].a = recompute_acceleration(n); // a(n)
		gas[n].r += gas[n].vr*dt + gas[n].a*(0.5*dt2); // x(n) -> x(n+1)

        if (gas[n].r < 0.)
          gas[n].r = 0.;

		const double a_old = gas[n].a;
		const double a_new = gas[n].a = recompute_acceleration(n); // a(n+1)
		gas[n].vr += (a_new+a_old)*0.5*dt; // v(n+1) = v(n) + 0.5*(a(n+1)+a(n))*dt

        gas[n].t += dt;
	}
}

void leapfrog_update()
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
		gas[n].a = recompute_acceleration(n); // a(n+1)
		gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 

        if (gas[n].r < 0.)
          gas[n].r = 0.;
        gas[n].t += dt;
	}

}

void adaptive_leapfrog_update()
{
    // Leapfrog (technically this should be equivalent to verlet)
    // v(n+1/2) = v(n) + a(n)*dt/2
    // x(x+1) = x(n) + v(n+1/2)*dt
    // v(n+1) = v(n+1/2) + a(n+1)*dt/2
    // Only the integer step velocity and positions are stored (not half step)
    // The time step is updated in a dumb way below - this version is not time-reversible/symplectic
    dt = get_time_step();
    dt2 = dt*dt;
    //cout << "New timestep = "<< dt << endl;
    #pragma omp parallel for
	for (size_t n = 0; n < gas.size(); ++n)
    {
		gas[n].vr += gas[n].a*dt/2; // v(n) -> v(n+1/2) 
		gas[n].r += gas[n].vr*dt; // x(n) -> x(n+1)
		gas[n].a = recompute_acceleration(n); // a(n+1)
		gas[n].vr += gas[n].a*dt/2; // v(n+1/2) -> v(n+1) 

        if (gas[n].r < 0.)
          gas[n].r = 0.;
        gas[n].t += dt;
	}

}
/*
inline void drift(double tau)
{
  for (size_t n = 0; n < gas.size(); ++n)
	  gas[n].r += gas[n].vr*tau; 
}

inline void kick(double tau)
{
	for (size_t n = 0; n < gas.size(); ++n)
    {
		gas[n].a = recompute_acceleration(n); 
		gas[n].vr += gas[n].a*tau;
    }
}

double select()
{
  double smallest_t_dyn;
  smallest_t_dyn = gas[0].t_dyn;
  for(size_t i = 0; i < nshells; ++i)
  {
    if (gas[i].t_dyn < smallest_t_dyn)
      smallest_t_dyn = gas[i].t_dyn;
  }
  return std::min(std::max(c_dyn*smallest_t_dyn, dt_min), dt_max);
}

void adaptive_reflexive_leapfrog_update(double tau_)
{
  // Recursive algorithm based on https://arxiv.org/pdf/astro-ph/9710043.pdf
  //
  double tau;
  tau = tau_;

  drift(tau);
  select();
  if (tau == dt_min)
  {
    kick(tau);
    drift(tau);
  }
  else{
    drift(-tau);
    adaptive_reflexive_leapfrog_update(tau/2.);
    kick(tau);
    adaptive_reflexive_leapfrog_update(tau/2.);
  }
}

*/
void output_shell_evolution()
{
  size_t i;
  cout << std::scientific; 
  double t = gas[0].t;

  cout << t << "    ";
//  for (size_t m = 0; m < r_avg.size(); m++)
//    cout << get_mass_interior_to_r(r_avg[m])/get_mass_interior_bg(t, r_avg[m]) << "    ";
  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    //cout << get_mass_interior(r_out[m]) << "    " << get_mass_interior_bg(t, gas[r_out[m]].r) << "    ";
    cout << gas[i].r << "    " << gas[i].vr << "    " << get_delta_interior(t, i) << "    " << get_mass_interior(i) << "    ";

  }
  double t_dyn;
  for (size_t m = 0; m < r_out.size(); m++)
  {
    i = r_out[m];
    t_dyn = gas[i].t_dyn;
    cout << t_dyn/dt << "    ";
  }

  /*
  for (size_t m = 0; m < gas.size(); m++)
    cout << gas[m].r << "    ";
  */
  cout << endl;
}

void animate(int)
{
    //adaptive_leapfrog_update();
    leapfrog_update();

#ifdef USE_OPEN_GL
	glutPostRedisplay();
	glutTimerFunc(delay,animate,0);
#endif // USE_OPEN_GL	

    if ((sim_time%100 == 0)) output_shell_evolution();

	sim_time++;		// Increment simulation time

	// Exit if sim_time_max exceeded
	if (sim_time*dt > sim_time_max) exit(0);
	//if (sim_time > 100000) exit(0);
}

int main(int argc, char **argv)
  //
{

#ifdef USE_OPEN_GL
	glutInit(&argc,argv);
	glutInitDisplayMode (GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize (640, 640);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("So Many Shells");
	glutDisplayFunc(display);
	glutTimerFunc(delay,animate,0);
#endif // USE_OPEN_GL	

	initialize_gas();

#ifdef USE_OPEN_GL
	glutMainLoop();
#endif // USE_OPEN_GL	

#ifndef USE_OPEN_GL
    while(true)
      animate(0);
#endif// USE_OPEN_GL
}
