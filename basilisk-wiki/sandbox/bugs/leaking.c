/**
# Leaking boundaries

This is based on [the shock test](/src/test/shock.c). */

#include "saint-venant.h"

int LEVEL = 5;

/**
We define a new boundary for the cylinder. */

bid cylinder;

int main() {
  size (5.);
  G = 9.81;
  origin (-L0/2., -L0/2.);
  //  gradient = zero;
  init_grid (1 << LEVEL);
  run();
}

/**
We impose height and velocity on the left boundary. */

#define H0 3.505271526
#define U0 6.29033769408481

event init (i = 0) {

  /**
  The geometry is defined by masking and the initial step function is
  imposed. 

  If this is replaced by "convex" boundaries there seems to be no
  problem. */

  mask (sq(x + 0.5) + sq(y) < sq(0.5) ? cylinder : none);
  
  foreach() {
    h[] = (x <= -1 ? H0 : 1.);
    u.x[] = (x <= -1 ? U0 : 0.);
  }
}

event logfile (i++; t <= 0.3) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

/**

The mass is not conserved as it should.

~~~gnuplot mass and time
plot 'log' u 1:5 w l
~~~
*/

