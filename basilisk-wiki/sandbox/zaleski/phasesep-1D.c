/**
# Phase separation in 1D

This tests the VdW EOS and spinodal decomposition */

#include "grid/multigrid1D.h"
#include "vdw.h"

/**
We start with initial conditions
etc... */

#define LEVEL 5

#define rho_c 1
#define R_g 1 
#define theta 0.95
#define p_c 1

#define MU 1e-3

double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*(8/(3-rhop) - 3*rhop/theta);
}
  
int main()
{
  periodic (right); 
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-3;
  run();
}

event init (i = 0)
{
  lambda=0.01;
  foreach()
    {
      mu[] = MU;
      rho[] = rho_c + 0.4*sin(2*pi*x/L0);
      q.x[] = q.y[] = 0.;
    }
   boundary(all);
}


event logfile (i += 1000) {
  stats s = statsf (rho);
  fprintf (stderr, "%g %d %g %g %g %.8f\n", t, i, dt, s.min, s.max, s.sum);
}

event outputfile (t <= 4; t += 0.01) {
  foreach()
    {
      fprintf (stdout, "%g %g %g\n", x, rho[],q.x[]);
 //     if (rho[] < 0.)
 //	{
 //	  fprintf(stderr,"Stop: negative density\n");
 //        return 1; // stops when rho becomes unphysical. 
 //	}
    }
  fprintf (stdout, "\n\n");
}

/**
~~~gnuplot  Evolution of the density and velocity
stats 'out' nooutput
set key top left
plot 'out' index (STATS_blocks-2)  w l t "final", 'out' index 0 w l t "initial"
~~~
*/
