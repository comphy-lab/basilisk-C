/**
# Merging of two vortices (centered RK solver)

This test is similar to [/src/test/stream.c]() but uses a
centered Navier--Stokes solver. */

#define BGHOSTS 2
#include "nsrk.h"

#define MAXLEVEL 8

int main() {
  origin (-0.5, -0.5);
  init_grid (1 << MAXLEVEL);
  run();
}

event init (t = 0) {
  scalar psi[], omega[];
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);
  double dd = 0.1;
  foreach() {
    omega[] = (exp(-(sq(x - dd) + sq(y))/(dd/10.)) +
	       exp(-(sq(x + dd) + sq(y))/(dd/10.)));
    psi[] = 0.;
  }
  boundary ({psi,omega});
  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});
}

event logfile (t <= 30; t += 0.5) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);
  fprintf (ferr, "%g %d %g %g %g %d\n", t, i, dt, s.sum, s.max, mgp.i);
}

event movie (t += 0.2; t <= 30) {
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, linear = true, file = "vort.mp4");
  foreach()
    omega[] = level;
  output_ppm (omega, spread = 2, file = "level.mp4");
}

#if TREE
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif

/**
## Results

![Vorticity](vortex/vort.mp4)

![Grid](vortex/level.mp4)

~~~gnuplot Vorticity statistics
set size square
set grid
set xlabel 'time'
plot 'log' u 1:4 w l lw 2 t 's.sum', '' u 1:5 w l lw 2 t 's.max'
~~~
*/
