/**
# Gravity wave

A similar test to the [capillary wave](capwave.c) but for a pure
gravity wave, using the [reduced gravity](/src/reduced.h) approach.

We use a constant-resolution grid, the Navier--Stokes solver for
two-phase flows and reduced gravity. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "prosperetti-gravity.h"

/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */

double se = 0; int ne = 0;

int main() {

  /**
  The domain is 2x2 to minimise finite-size effects. The viscosity is
  constant. The acceleration of gravity is 50. */

  L0 = 2.;
  Y0 = -L0/2.;
  G.y = 50.;
  rho1 = 1, rho2 = 0.1;
  mu1 = mu2 = 0.0182571749236;
  TOLERANCE = 1e-6[*];

  /**
  We vary the resolution to check for convergence. */

  for (N = 16; N <= 128; N *= 2) {
    se = ne = 0;
    run();
  }
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */

event init (t = 0) {
  fraction (f, y - 0.01*cos (2.*pi*x));
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);

/**
We output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 0.00225584983639310905; t <= 1.66481717925811447992) {

  /**
  To get an accurate amplitude, we reconstruct interface position
  (using height functions) and take the corresponding maximum. */

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  
  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*16.032448313657, max);
  fflush (fp);

  /**
  To compute the RMS error, we get data from the reference file
  [prosperetti-gravity.h]() and add the difference to the accumulated
  error. */

  se += sq(max - prosperetti[ne][1]); ne++;

  if (N == 64)
    output_facets (f, stdout);
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0, sqrt(se/ne)/0.01);

/**
## Results

~~~gnuplot Evolution of the amplitude of the gravity wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../prosperetti-gravity.h' u 2:4 w l t "Prosperetti", \
     'wave-128' every 10 w p t "Basilisk"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
set grid
plot [5:200][1e-4:1]'log' u 1:2 t "Basilisk" w lp, 2./x**2 t "Second order"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html#gravity)
*/
