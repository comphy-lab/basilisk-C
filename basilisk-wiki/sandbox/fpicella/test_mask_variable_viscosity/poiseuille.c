/**
# Poiseuille flow in a periodic channel
but oriented in x direction.
I try to work with simple maks function
*/

#include "navier-stokes/centered.h"
#include "view.h"

double HEIGHT = 1.0;

int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4. [0]);
  DT = HUGE [0];
  
  origin (-L0/2., -L0/2.);
  periodic (right);
  //periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
  //for (N = 16; N <= 64; N *= 2)
	N=64;
    run();
}

scalar un[];

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {

  /**
  The gravity vector is aligned with the channel and viscosity is
  unity. */
  
  a[] = {1.,1.};
  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  
  mask (y > +HEIGHT/2. ? top : none);
  mask (y < -HEIGHT/2. ? bottom : none);
	
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);

  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.y[];
}

/**
We check for a stationary solution. */

event logfile (t += 0.1; i <= 1000) {
//  double du = change (u.y, un);
//  if (i > 0 && du < 1e-6)
//    return 1; /* stop */
}
