/**
# single rigid particle in a viscous fluid.
using Smoothed Profile method
*/

#define EPS 1e-8
double RADIUS = 0.25;
double force_tune = 0.5;
double ZITA       = 1.0;

#include "embed.h"
#include "navier-stokes/centered.h"
#include "fpicella/src/periodic-shift-treatment.h"
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/SmoothedProfile.h"

Particles colloid; // colloid is the name of the variable containing the particles I'm working with

#define DOMAIN -sq(x)-sq(y)+sq(0.5) 
//#define DOMAIN intersection(-y+0.5+EPS,y+0.5+EPS)


int main()
{
	display_control(force_tune,-100,100);
	display_control(ZITA,0.001, 100);

  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (2. [0]);
  DT = 1.0 [0];
  
  origin (-L0/2., -L0/2.);

	/**
	Set variable for variable forces */
	a = av;
  
  stokes = true;
	periodic(left);
	periodic(top);
	N=64;
  run();
}

event init (t = 0) {

  /**
  The viscosity is unity. */
  
  mu = fm;
	
  /**
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */

  solid (cs, fs, DOMAIN);

  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
  
  u.n[embed] = dirichlet (-y);
  u.t[embed] = dirichlet (+x);

}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000000) {

}

/**
Initialize colloid location. */
event init(i=0)
{
	colloid = init_tp_circle(1);
	foreach_particle(){
		p().x = 0.+L0/2.*0.0;
		p().y = 0.0;
		fprintf(stderr,"toto %+6.5e %+6.5e %d \n",p().x,p().y,_l_particle);
		p().u.x = 0.;
		p().u.y =-0.;
		p().w.x = -0.0;
		p().w.y = -0.0;
	}
}
