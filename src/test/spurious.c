/**
# Circular droplet in equilibrium

This is the classical "spurious" or "parasitic currents" test case
discussed in [Popinet, 2009](/src/references.bib#popinet2009).

We use the Navier--Stokes solver with VOF interface tracking and
surface tension. */

#define JACOBI 1

#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

/**
The interface is represented by the volume fraction field *c*. */

scalar c[], * interfaces = {c};

/**
The diameter of the droplet is 0.8. The density is constant (equal to
unity by default), and the viscosity is defined through the Laplace
number
$$
La = \sigma\rho D/\mu^2
$$
with $\sigma$ set to one. The simulation time is set to the
characteristic viscous damping timescale. */

#define DIAMETER 0.8
#define MU sqrt(DIAMETER/LAPLACE)
#define TMAX (sq(DIAMETER)/MU)

/**
We will vary the number of levels of refinement (to study the
convergence), the Laplace number and *DC* a convergence parameter
which measures the variation in volume fraction between successive
timesteps (to evaluate whether we are close to a steady solution). */

int LEVEL;
double LAPLACE;
double DC = 0.;
FILE * fp = NULL;

int main() {
  
  /**
  We neglect the advection terms and vary the Laplace, for a constant
  resolution of 5 levels. */

  DT = HUGE [0];
  TOLERANCE = 1e-6 [*];
  stokes = true;
  c.sigma = 1;
  LEVEL = 5;
  N = 1 << LEVEL;
  for (LAPLACE = 120; LAPLACE <= 12000; LAPLACE *= 10)
    run();

  /**
  We now fix the Laplace number and look for stationary solutions
  (i.e. the volume fraction field is converged to within 1e-10) for
  varying spatial resolutions. */

  LAPLACE = 12000; DC = 1e-10;
  for (LEVEL = 3; LEVEL <= 7; LEVEL++) 
    if (LEVEL != 5) {
      N = 1 << LEVEL;
      run();
    }
}

/**
We allocate a field to store the previous volume fraction field (to
check for stationary solutions). */

scalar cn[];

event init (i = 0) {

  /**
  We set the constant viscosity field... */

  mu[] = {MU,MU};

  /**
  ... open a new file to store the evolution of the amplitude of
  spurious currents for the various LAPLACE, LEVEL combinations... */

  char name[80];
  sprintf (name, "La-%g-%d", LAPLACE, LEVEL);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");

  /**
  ... and initialise the shape of the interface and the initial volume
  fraction field. */
  
  fraction (c, sq(DIAMETER/2) - sq(x) - sq(y));
  foreach()
    cn[] = c[];
}

event logfile (i++; t <= TMAX)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  
  double dc = change (c, cn);
  if (i > 1 && dc < DC)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g %g\n",
	   MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER), dc);
}

event error (t = end) {
  
  /**
  At the end of the simulation, we compute the equivalent radius of
  the droplet. */

  double vol = statsf(c).sum;
  double radius = sqrt(4.*vol/pi);

  /**
  We recompute the reference solution. */
  
  scalar cref[];
  fraction (cref, sq(DIAMETER/2) - sq(x) - sq(y));
  
  /**
  And compute the maximum error on the curvature *ekmax*, the norm of
  the velocity *un* and the shape error *ec*. */
  
  double ekmax = 0.;
  scalar un[], ec[], kappa[];
  curvature (c, kappa);
  foreach() {
    un[] = norm(u);
    ec[] = c[] - cref[];
    if (kappa[] != nodata) {
      double ek = fabs (kappa[] - (/*AXI*/ + 1.)/radius);
      if (ek > ekmax)
	ekmax = ek;
    }
  }
  
  /**
  We output these on standard error (i.e. the *log* file). */

  norm ne = normf (ec);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", 
	   LEVEL, LAPLACE, 
	   normf(un).max*sqrt(DIAMETER), 
	   ne.avg, ne.rms, ne.max,
	   ekmax);
}

/**
We use an adaptive mesh with a constant (maximum) resolution along the
interface. */

#if TREE
event adapt (i <= 10; i++) {
  adapt_wavelet ({c}, (double[]){0}, maxlevel = LEVEL, minlevel = 0);
}
#endif

/**
## Results

The maximum velocity converges toward machine zero for a wide range of
Laplace numbers on a timescale comparable to the viscous dissipation
timescale, as expected.

~~~gnuplot Evolution of the amplitude of the capillary currents $\max(|\mathbf{u}|)(D/\sigma)^{1/2}$ as a function of non-dimensional time $\tau=t\mu/D^2$ for the range of Laplace numbers indicated in the legend.
set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", \
  'La-12000-5' w l t "La=12000"
~~~

The equilibrium shape and curvature converge toward the exact shape
and curvature at close to second-order rate.

~~~gnuplot Convergence of the error on the equilibrium shape of the droplet with resolution. The diameter is given in number of grid points.
set xlabel 'D'
set ylabel 'Shape error'
set logscale x
set xtics 2
set pointsize 1
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):5 w lp t "RMS", \
            '< sort -n -k1,2 log' u (0.8*2**$1):6 w lp t "Max", \
             0.2/(x*x) t "Second order"
~~~

~~~gnuplot Convergence of the relative error on the equilibrium curvature value with resolution. The diameter is given in number of grid points.
set ylabel 'Relative curvature error'
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):($7/2.5) w lp t "Max", \
             0.6/(x*x) t "Second order"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/
