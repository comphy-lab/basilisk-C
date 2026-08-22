/**
# Circular droplet in equilibrium

This is the classical "spurious" or "parasitic currents" test case
discussed in [Popinet, 2009](/src/references.bib#popinet2009).

We use the Navier--Stokes solver with VOF interface tracking and
surface tension. */

#define GAUSS_SEIDEL 1
#define LB_HP_FILTER 0

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension-ssvm.h"

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
#define TMAX max(1.2*sq(DIAMETER)/MU, 100)

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

  size(1);
  //origin(-L0/2, -L0/2);
  double ratios[2] = {1., 15.};
  ratio_dt_st = 15.;

  DT = HUGE [0];
  TOLERANCE = 1e-10 [*];
  //TOLERANCE_MU = 1e-6;
  f.sigma = 1;
  stokes = true;
  LEVEL = 6;
  LAPLACE = 12000; DC = 1e-10;
  mu1 = mu2 = MU;
  N = 1 << LEVEL;

  for (int i = 0; i <= 1; i++) {
    ratio_dt_st = ratios[i];
    run();
  }

  return 1;
}

scalar kappa_test2[];

event end_timestep(i++){
 fprintf(fout, "stability %12e %12e %12e\n", dt, sqrt(1.0*pow(L0/(1 << grid->maxdepth), 3)/M_PI/f.sigma), MU/f.sigma*L0/(1 << grid->maxdepth));
 //weight_kappa = (double)(i%10 == 0 ? 1 : 0.2);
 curvature (f, kappa_test2, 1, add = false);
}

/**
We allocate a field to store the previous volume fraction field (to
check for stationary solutions). */

scalar cn[];

event init (i = 0) {

  //mu[] = {MU,MU};

  /**
  ... open a new file to store the evolution of the amplitude of
  spurious currents for the various LAPLACE, LEVEL combinations... */

  char name[80];
  sprintf (name, "La-%g-%d-%d", LAPLACE, LEVEL, (int) ratio_dt_st);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");

  /**
  ... and initialise the shape of the interface and the initial volume
  fraction field. */
  
  fraction (f, DIAMETER/2 - sqrt(sq(x) + sq(y)));
  foreach()
    cn[] = f[];
}

event logfile (i++; t <= TMAX)
// event logfile (i++; i <= 200)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  
  double dc = change (f, cn);
  //if (i > 1 && dc < DC)
  //  return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g %g %g\n",
	   MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER), normf(un).rms*sqrt(DIAMETER), dc);
  fflush(fp);
}

event error (t = end) {
  
  /**
  At the end of the simulation, we compute the equivalent radius of
  the droplet. */

  double vol = statsf(f).sum;
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
  curvature (f, kappa);
  foreach(reduction(max:ekmax)) {
    un[] = norm(u);
    ec[] = f[] - cref[];
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

/*
#if TREE
event adapt (i <= 10; i++) {
  adapt_wavelet ({c}, {0}, maxlevel = LEVEL, minlevel = 0);
}
#endif
*/

/**
## Results

The maximum velocity converges toward machine zero for a wide range of
Laplace numbers on a timescale comparable to the viscous dissipation
timescale, as expected.

Evolution of the amplitude of the capillary currents $\max(|\mathbf{u}|)(D/\sigma)^{1/2}$ as a function of non-dimensional time $\tau=t\mu/D^2$ for the range of timesteps.

~~~gnuplot Evolution of the amplitude of the capillary currents
set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
plot 'La-12000-6-1' w l t "{/Symbol D} t = {/Symbol D} t_{ST}, CSVM", \
  'La-12000-6-15' w l t "{/Symbol D} t = 15 {/Symbol D} t_{ST}, CSVM"
~~~

## See also

* [Same test with Gerris](https://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/
