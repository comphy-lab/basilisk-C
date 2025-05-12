/**
# Bouncing Saint-Venant bump with the Mc Cormack scheme in 1D.

This test case is a restriction to 1D of [bump2D-vdw.c]() for debugging purposes */

#include "grid/cartesian1D.h"
#include "vdw.h"

/**
We start with initial conditions
etc... as when using the standard 
[bump2D.c](/src/test/bump2D.c) . */

#define LEVEL 8

double P0(double x)
{
  return 0.5*x*x;
}
  
int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-4;
  run();
}

event init (i = 0)
{
  lambda=0.;
  foreach()
    {
      mu[] = 1.e-3;
      rho[] = 0.1 + exp(-200.*(x*x + y*y)/(L0*L0));
      q.x[] = q.y[] = 0.;
    }
}


event logfile (i++) {
  stats s = statsf (rho);
  fprintf (stderr, "%g %d %g %g %g %.8f\n", t, i, dt, s.min, s.max, s.sum);
}

event outputfile (t <= 0.4; t += 0.2) {
  foreach()
    {
      fprintf (stdout, "%g %g %g\n", x, rho[],q.x[]);
      //if (rho[] < 0.)
         //return 1; // stops when rho becomes unphysical. 
    }
  fprintf (stderr, "\n");
}

/**

# Same test with Saint-Venant solver 
 see [bump-saint-venant1D.c]()


# Results
~~~gnuplot  Evolution of the density and momentum. The results are comparable to that obtained when running with the Saint-Venant solver but the MacCormack scheme has oscillations.
plot 'out' u 1:2 w l t "rho", 'out' u 1:3 w l t "q.x"
~~~

*/
