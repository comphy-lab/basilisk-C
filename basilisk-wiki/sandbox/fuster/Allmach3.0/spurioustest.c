/**
# Circular droplet in equilibrium

This is the classical "spurious" or "parasitic currents" test case
discussed in [Popinet, 2009](/src/references.bib#popinet2009)
to test the compressible Navier--Stokes solver with VOF interface tracking and
surface tension. */

#define JACOBI 1

#include "two-phase-compressible.h"
#include "compressible-tension.h"

#define DIAMETER 0.8
#define MU sqrt(DIAMETER/LAPLACE)
/*Change 0.1 to 5 to obtain the results plotted at the end of this page */
#define TMAX (5.*sq(DIAMETER)/MU)

int LEVEL;
double LAPLACE;
double DC = 0.;
FILE * fp = NULL;

int main() {
  
  /**
  We neglect the advection terms and vary the Laplace, for a constant
  resolution of 5 levels. */
  
  PI1 = 300.;
  gamma1 = 7.14;
  gamma2 = 1.4;

  TOLERANCE = 1e-6 [*];
  L0 = 1. [0];
  DT = HUGE [0];
  
  f.sigma = 1;
  f.gradient = zero;
  LEVEL = 5;
  N = 1 << LEVEL;
  for (LAPLACE = 120; LAPLACE <= 12000; LAPLACE *= 10) {

    /**
    We set the constant viscosity field... */

    mu1 = mu2 = MU;
    run();
  }
}

/**
We allocate a field to store the previous volume fraction field (to
check for stationary solutions). */

scalar cn[];

event init (i = 0) {

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
  
  double p0L = 1.;
  double p0 = p0L + f.sigma/DIAMETER*2;
  fraction (f, sq(DIAMETER/2) - sq(x) - sq(y));
  foreach() {
    cn[] = f[];
    frho1[]  = f[];
    frho2[]  = (1. - f[]);
    double pL = p0L;
    p[]   = pL*f[] + p0*(1.-f[]);
    fE1[]   = f[]*(pL/(gamma1 - 1.) + PI1*gamma1/(gamma1 - 1.));
    fE2[]   = (1.-f[])*p0/(gamma2 - 1.);
    q.x[] = 0.;
    q.y[] = 0.;
  }
  boundary ({cn});
}

event logfile (i++; t <= TMAX)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  
  double dc = change (f, cn);
  if (i > 1 && dc < DC)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach()
    un[] = norm(q);
  fprintf (fp, "%g %g %g\n",
     MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER), dc);
}

/**
## Results

The maximum velocity converges toward machine zero for a wide range of
Laplace numbers on a timescale comparable to the viscous dissipation
timescale, as expected.

~~~gnuplot Evolution of the amplitude of the capillary currents $\max(|\mathbf{u}|)(D/\sigma)^{1/2}$ as a function of non-dimensional time $\tau=t\mu/D^2$ for the range of Laplace numbers indicated in the legend.
set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", 'La-12000-5' w l t "La=12000"
~~~

The high frequency oscillations scale with the bubble resonance frequency $\omega_0 = \frac{1}{R_0}\sqrt{\frac{3 \gamma p_0}{\rho_l}}$

~~~gnuplot  In the unit system used this implies $\tau'=t \omega_0/(2 \pi)$ 
set xlabel 't f_0'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
DIAMETER=0.8
plot [0:10] 'La-120-5' u ($1*sqrt(120*3*1.4*(1. + 2./DIAMETER))/(pi*DIAMETER)):2 w l t "La=120", \
     'La-1200-5' u ($1*sqrt(1200*3*1.4*(1. + 2./DIAMETER))/(pi*DIAMETER)):2 w l t "La=1200", \
     'La-12000-5' u ($1*sqrt(12000*3*1.4*(1 + 2./DIAMETER))/(pi*DIAMETER)):2 w l t "La=12000"
~~~

## See also

* [Convergence study of the same test](http://basilisk.fr/sandbox/fuster/Allmach3.0/spurioustest2.c)
* [Same test with the incompressible version](http://basilisk.fr/src/test/spurious.c)
* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/
