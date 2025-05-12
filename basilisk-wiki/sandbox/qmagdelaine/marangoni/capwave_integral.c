/**
# Capillary wave

This is the classical test case first proposed in [Popinet & Zaleski,
1999](/src/references.bib#popinet1999).

We use a constant-resolution grid, the Navier--Stokes solver with VOF
interface tracking and surface tension. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"

/**
We use the integral formulation of the surface tension. */

#include "integral_vof_quasi1D.h"
#include "test/prosperetti.h"

/**
The interface is represented by the volume fraction field *c*. */

scalar c[], * interfaces = {c};

/**
In the integral formulation the surface tension can be a variable, but here
we set to 1 everywhere and we declare liquid_sigma[] to store it. */

void surface_tension (scalar gamma) {
  foreach()
    gamma[] = 1.;
  boundary({gamma});
}
scalar liquid_sigma[];

/**
As with [contact.h](/src/contact.h), we define vector field to store the height
functions assciated to the VOF tracer f. */

vector h[];

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

  L0 = 2.;
  Y0 = -L0/2.;
  TOLERANCE = 1e-6;
  const face vector muc[] = {0.0182571749236, 0.0182571749236};
  mu = muc;

  /**
  We associate the height functions $h$ to the VOF tracer $c$. */

  c.height = h;

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

  /**
  With the current version of [intergral_vof_quasi1D.h](intergral_vof_quasi1D.h),
  the phase corresponding to $f = 1$ has to be the lower phase. It was not the
  case before, so it is probably just a typo in the header file. It will be
  corrected soon. */

  fraction (c, 0.01*cos (2.*pi*x) - y);
  
  /**
  We set the surface tension. */
  
  surface_tension (liquid_sigma);
  c.sigma = liquid_sigma;
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
  position (c, pos, {0,1});
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
set output 'amplitude.png'
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../prosperetti.h' u 2:4 w l t "Prosperetti", \
     'wave-128' every 10 w p t "Basilisk"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set output 'convergence.png'
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
set grid
plot [5:200][1e-4:1]'log' t "Basilisk" w lp, 2./x**2 t "Second order"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html)
*/