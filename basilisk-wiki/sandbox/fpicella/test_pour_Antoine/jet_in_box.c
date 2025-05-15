/**
# Poiseuille flow in a periodic channel
but oriented in x direction.
I try to work with simple maks function
*/

#include "navier-stokes/centered.h"
#include "view.h"

double HEIGHT = 100.;
double b = 1.0;
double epsilon = 0.1;

double Reynolds = 250.;

face vector muv[]; // will contain the problem's viscosity.

int main()
{

	display_control(Reynolds,0.01,1000);
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (16.);
  DT = 1.0;
  
  origin (-L0/2., -L0/2.);
  //periodic (right);
  //periodic (top);
  
  //stokes = true;
  //TOLERANCE = 1e-7;
  
  //for (N = 16; N <= 64; N *= 2)
	mu = muv;
	
	N=128;
    run();
}

/** 
	Define the boundary associated to the _plate. */
bid plate;

event init (t = 0) {

  /**
  The gravity vector is aligned with the channel and viscosity is
  unity. */
  
//  a[] = {1.,1.};
//  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  
  mask (y > +HEIGHT/2. ? top : none);
  mask (y < -HEIGHT/2. ? bottom : none);
	/**
	Define the location of the _plate_, where the jet will impinge */
	//mask ((x > 0. && x <L0/2.-epsilon*10.) && fabs(y) < epsilon ? plate : none);
	mask ((x > 0.) && fabs(y) < epsilon ? plate : none);
	
	u.n[plate] = dirichlet(0.);
	u.t[plate] = dirichlet(0.);
	
	u.n[left] = dirichlet( fabs(y)>b/2. ? 0 : 1.); 
	p[left]    = neumann(0.);
	pf[left]   = neumann(0.);

	u.n[right] = dirichlet(1./L0);
	p[right]   = neumann(0.);
	pf[right]  = neumann(0.);

	//u.n[right] = neumann(0.);
	//p[right]   = dirichlet(0.);
	//pf[right]  = dirichlet(0.);

}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
}

/**
We check for a stationary solution. */
scalar omega[];

event logfile (t += 1.0; i <= 100000) {
	vorticity(u,omega);
//  double du = change (u.y, un);
//  if (i > 0 && du < 1e-6)
//    return 1; /* stop */
}
