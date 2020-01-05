#pragma once

#include <iostream>
using std::ostream;

#include <iomanip>
using std::setw;

#include <cmath>
using std::sqrt;

class shell
{
public:
	// member data
	double mass;
	double r;  // position
	double vr; // radial velocity
	double a;  // radial acceleration
    double l;  // angular momentum
	double t;  // current time;
	double t_dyn;  // characteristic dynamical time;

	// constructor
	shell(double m_, double r_, double vr_, double l_) : mass(m_), r(r_), vr(vr_), a(0.), l(l_), t(0.) , t_dyn(0.) {}
	shell(double m_, double r_, double vr_, double l_, double t_) : mass(m_), r(r_), vr(vr_), a(0.), l(l_), t(t_), t_dyn(0.) {}
	
	// methods (or member functions)
	void evolve(void (&method)(shell&,double), double dt)
	{
		// *this refers to the active shell object
		method(*this, dt);
	}

#ifdef USE_OPEN_GL
	void draw(void)
	{
        const float radius = r;
        glBegin(GL_LINE_LOOP); 
          glColor3f(1.0, 0.5, 0.25); 
          for(int j = 0; j < 360; ++j)
            glVertex2f(radius*cos(j*M_PI/180.0f),radius*sin(j*M_PI/180.0f));
        glEnd();
		glFlush();
	}
#endif // USE_OPEN_GL	
	
	// a friend function does not belong to the class but is affiliated with it
	friend ostream& operator<<(ostream& os, const shell &s)
	{
		os.precision(8);
		os << setw(20) << s.t
		    << setw(20) << s.r
		    << setw(20) << s.vr
		    << setw(20) << s.a;
		return os;
	}
};
