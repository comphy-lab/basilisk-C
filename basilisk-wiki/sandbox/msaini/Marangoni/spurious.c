/**
# Circular droplet in equilibrium

This is the classical "spurious" or "parasitic currents" test case
discussed in [Popinet, 2009](/src/references.bib#popinet2009). We will to check the well-balancedness of the new [CSS method](/src/integral.h) and compare with the standard [CSF method](/src/tension.h) available in basilisk. */

#define JACOBI 1
//#define CSF 1


#include "navier-stokes/centered.h"
#if CSF == 1
#include "vof.h"
#include "tension.h"
scalar f[], * interfaces = {f};
#else
#include "two-phase-clsvof.h"
#include "integral.h"
scalar sigmaf[];
#endif

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
We will vary the Laplace number and check the spurious currents in the domain. */

int LEVEL;
double LAPLACE = 120.;
double DC = 0.;
FILE * fp = NULL;

int main() {
  
  /**
  We neglect the advection terms and vary the Laplace, for a constant
  resolution of 6 levels. */

  TOLERANCE = 1e-6 [*];
  stokes = true;
#if CSF == 1
  f.sigma = 1;
#else
  d.sigmaf = sigmaf;
#endif
  
  LEVEL = 6;
  N = 1 << LEVEL;
  for (LAPLACE = 1200; LAPLACE <= 12000; LAPLACE *= 10)
    run();
}

event init (i = 0) {

  /**
  We set the constant viscosity field... */

  const face vector muc[] = {MU,MU};
  mu = muc;

  /**
  ... open a new file to store the evolution of the amplitude of
  spurious currents for the various LAPLACE... */

  char name[80];
#if CSF == 1
  sprintf (name, "La-%g", LAPLACE);
#else
  sprintf (name, "La-CSS-%g", LAPLACE);
#endif
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");

  /**
  ... and initialise the shape of the interface and the initial volume
  fraction field. */
  
#if CSF == 1
  fraction (f, sq(DIAMETER/2) - sq(x) - sq(y));
#else
  foreach(){
    d[] = - sqrt (sq(x) + sq(y)) + DIAMETER/2.;
    sigmaf[] = 1.;
  }
#endif
}

event logfile (i++; t <= TMAX)
{
  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g\n",MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER));
}

event error (t = end) {
  
  scalar un[];
  foreach() {
    un[] = norm(u);
    }
  
  fprintf (stderr, "%d %10.10f %10.10f\n", 
           LEVEL, LAPLACE,normf(un).max*sqrt(DIAMETER));
}


/**
We use an adaptive mesh with a constant (maximum) resolution along the
interface. */

#if TREE
event adapt (i <= 10; i++) {
  adapt_wavelet ({f}, (double[]){0}, maxlevel = LEVEL, minlevel = 0);
}
#endif

/**
## Results

In the case of CSF model, the maximum velocity decays toward machine zero for a wide range of Laplace numbers on a timescale comparable to the viscous dissipation timescale, as expected. However, for the new CSS model, the velocity does not decay.

~~~gnuplot Comparison of two methods

set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
set key center right font "times,12"

plot 'La-120' w l t "CSF,La=120",\
     'La-1200' w l t "CSF,La=1200",\
     'La-12000' w l t "CSF,La=12000",\
     'La-CSS-120' w l t "CSS,La=120",\
     'La-CSS-1200' w l t "CSS,La=1200",\
     'La-CSS-12000' w l t "CSS,La=12000"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/