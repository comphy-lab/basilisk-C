/**
# An `ascii`-art example

Have a look at [the output stream](taa/out). 

This move is generated for addional context:

![A dipole-wall collision](taa/vor.mp4)

*/
#include "navier-stokes/centered.h"
#include "ascii-art.h"

double Re = 2500;

u.t[bottom] = dirichlet (0);

int main() {
  L0 = 15;
  X0 = Y0 = -L0/2;
  const face vector muc[] = {1./Re, 1./Re};
  mu = muc;
  N = 128;
  run();
}

event init (t = 0) {
  foreach() 
    u.y[]  = -exp(-sq(x) - sq(y));
  boundary ({u.y});
}

event art (i += 5) {
  scalar vor[];
  vorticity (u, vor);
  foreach()
    vor[] = fabs(vor[]);
  ascii_art(vor, min = 0, max = 1, n = 60);
}

event mov (i += 5) {
  scalar vor[];
  vorticity (u, vor);
  foreach()
    vor[] = fabs(vor[]);
  output_ppm (vor, file = "vor.mp4", n = 500,
	      min = 0, max = 1, map = gray);
}

event adapt (i++) 
  adapt_wavelet({u.x, u.y}, (double[]){0.002, 0.002}, 9);

event stop (t = 100);
