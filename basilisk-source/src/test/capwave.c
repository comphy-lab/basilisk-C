/**
# Capillary wave

This is the classical test case first proposed in [Popinet & Zaleski,
1999](/src/references.bib#popinet1999).

We use a constant-resolution grid, the Navier--Stokes solver with VOF
interface tracking (optionally coupled with levelset) and surface
tension (optionally using the integral formulation). */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#if CLSVOF
# include "two-phase-clsvof.h"
# include "integral.h"
# include "curvature.h"
#else
# include "vof.h"
# include "tension.h"
/**
The interface is represented by the volume fraction field *c*. */

scalar f[], * interfaces = {f};
#endif
#include "prosperetti.h"

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
  The domain is 2x2 to minimise finite-size effects. The surface
  tension is one and the viscosity is constant. */

  size (2. [1]);
  Y0 = -L0/2.;
#if CLSVOF
  const scalar sigma[] = 1.;
  d.sigmaf = sigma;
#else
  f.sigma = 1.;
#endif
  TOLERANCE = 1e-6 [*];
  mu[] = {0.0182571749236, 0.0182571749236};

  
  /**
  We vary the resolution to check for convergence. */

  for (N = 16; N <= 128; N *= 2) {
    se = 0, ne = 0;
    run();
  }
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */

event init (t = 0) {
  double k = 2., a = 0.01;
#if CLSVOF
  foreach()
    d[] = y - a*cos (k*pi*x);
#else
  fraction (f, y - a*cos (k*pi*x));
#endif
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);

/**
We output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 3.04290519077e-3; t <= 2.2426211256) {

  /**
  To get an accurate amplitude, we reconstruct interface position
  (using height functions) and take the corresponding maximum. */

  scalar pos[];
  position (f, pos, {0,1 [0]});
  double max = statsf(pos).max;

  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*11.1366559937, max);
  fflush (fp);

  /**
  To compute the RMS error, we get data from the reference file
  *prosperetti.h* and add the difference to the accumulated error. */

  se += sq(max - prosperetti[ne][1]); ne++;
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0, sqrt(se/ne)/0.01);

/**
## Results

~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../prosperetti.h' u 2:4 w l t "Prosperetti", \
     'wave-128' every 10 w p t "Basilisk", \
     '../capwave-clsvof/wave-128' every 10 w p t "Basilisk (CLSVOF)"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
set grid
plot [5:200][1e-4:1]\
     'log' t "Basilisk" w lp, 2./x**2 t "Second order", \
     '../capwave-clsvof/log' t "Basilisk (CLSVOF)" w lp, \
     2./x**2 t "Second order"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html)
*/
