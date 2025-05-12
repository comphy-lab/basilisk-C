/**
# Capillary wave for fluids with different densities

We use a constant-resolution grid, the Navier--Stokes solver with VOF
interface tracking and surface tension. */
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "prosperetti-density.h"


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
  f.sigma = 1;
  TOLERANCE = 1e-6;
  mu1= 0.0182571749236;
  mu2= 0.0182571749236;
  /**
  the density of fluid 1 is one and the density ratio is 10 */
  rho1=1;
  rho2=0.1;
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
  To get an accurate amplitude, we reconstruct the height function
  field and take the corresponding maximum. */

  vector h[];
  heights (f, h);
  double max = - HUGE;;
  foreach() 
    if (f[] > 0 && f[] < 1) {
      double yi = y + height(h.y[])*Delta;
      if (yi > max)
	max = yi;
    }

  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*15.016663878457, max);
  fflush (fp);

  /**
  To compute the RMS error, we get data from the reference file
  *prosperetti-density.h* and add the difference to the accumulated error. */

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
plot '../prosperetti-density.h' u 2:4 w l t "Prosperetti", \
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

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html#density)
*/