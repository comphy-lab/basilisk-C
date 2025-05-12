/**
# (Naive) Forced turbulence in a 2D-periodic box with non-coalescing droplets

We investigate the dynamics of droplets in a two-dimensional periodic box. To induce turbulence, we employ a velocity linearly dependent forcing term, similar to the approach described in http://basilisk.fr/src/examples/isotropic.c. It should be noted that the intention is not to specifically generate two-dimensional isotropic turbulence with this forcing term, but rather to provide a crude way to make the droplets collide. To prevent numerical coalescence of droplets, we utilize the no-coalescence module. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "no-coalescence.h"


/**
We need to store the variable forcing term. */

face vector av[];


int maxlevel = 8;

int main (int argc, char * argv[])
{
  //maxruntime (&argc, argv);
  if (argc > 1)
    maxlevel = atoi(argv[1]);

  /**
  The domain is $(2\pi)^3$ and doubly-periodic. */
  
  L0 = 2.*pi;
  foreach_dimension()
    periodic (right);

  /**
  The acceleration is defined below. The
  level of refinement is *maxlevel*. */

  double rhoInt = 1.;
  double muInt = 0.01;

  rho1 = rhoInt, mu1 = muInt;
  rho2 = rhoInt, mu2 = muInt;
  f.sigma = 4.;
  
  a = av;
  N = 1 << maxlevel;
  run();
}

/**
## Initial positions of the drop */

double geometry (double x, double y) {

  double inout=0.,i1=0.,i2=0.,i3=0.,i4=0.,Circle1=0.,Circle2=0.,Circle3=0.,Circle4=0.;

  Circle1 = -sq(x-pi/2)-sq(y-pi/2)+sq(1);
  Circle2 = -sq(x-pi/2)-sq(y-3.5*pi/2)+sq(1);
  Circle3 = -sq(x-3.5*pi/2)-sq(y-3*pi/2)+sq(1);
  Circle4 = -sq(x-3*pi/2)-sq(y-pi/2)+sq(1);

  if(Circle1 > 0.) i1=1;
  if(Circle2 > 0.) i2=1;
  if(Circle3 > 0.) i3=1;
  if(Circle4 > 0.) i4=1;

  return(i1+i2+i3+i4);
}


/**
## Initial conditions

The initial condition is a periodic flow with white noise. */

event init (i = 0) {
  if (!restore (file = "restart"))

    fraction(f, geometry(x,y));

    foreach() {
      u.x[] = (1-f[])*cos(pi*y) + noise()*0.1;
      u.y[] = (1-f[])*sin(pi*x) + noise()*0.1;
    }

}

/**
## Linear forcing term

We compute the average velocity and add the corresponding linear
forcing term. */

event acceleration (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  foreach_face()
    av.x[] += 0.1*(1-f[])*((u.x[] + u.x[-1])/2. - ubar.x);
}


event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
  
}

/**
We generate movies of the vortices, volume fraction and velocity (along x). 
*/

event output (t += 0.025; t <= 15.)
{
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, min = -1, max = 1, file = "omega.mp4");
  output_ppm (f, file = "fraction.mp4");
  output_ppm (u.x, file = "u.mp4");
}


