/**
# Rising Gaussian plume test

This rising Gaussian plume test consists of multiple files and scripts so that we can verify the adaptive grid convergence using a reference solution obtained with a high-resolution equidistant grid. 

1. [rise.c]() is used to generate the reference-solution dump files
2. [risea.c]() runs the case with an adaptive grid
3. [genarr_rise.c]() converts the equidistant dumps to 2D arrays of the solution at levels
4. [get-error-rise.c]() Compares the adaptive solution to the reference arrays
5. [drawit.c]() generates plots of the solution
*/

#include "nsf4t.h"
scalar b[], * tracers = {b};
scalar zeros[]; // Just zeros

int main(int argc, char ** argv) {
  foreach_dimension()
    periodic (left);
  int l = 9;
  if (argc >= 2) 
    l = atoi(argv[1]);
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2;
  const scalar muz[] = 0.001;
  kappa = nu = muz;
  N = 1 << l;
  a.x = zeros;
  a.y = b;
  run();
}

event init (t = 0) {
  TOLERANCE = 1e-10;
  foreach_vert() 
    b[] = exp (-sq((x - exp(1)/10.)) - sq(y + 1) - sq(z - pi/10.)) ;
  boundary ({b});
}

event mov (t += 0.1) 
  output_ppm (b, file = "b.mp4", n = 300);

event dumper (t = {1,3,5}) {
  char str[99];
  restriction({b});
  sprintf (str, "%d-dump%d", (int)(t + 0.1), depth());
  dump(str, list = {b});
}

