/**
# Circular droplet in equilibrium

This is the classical "spurious" or "parasitic currents" test case
discussed in [Popinet, 2009](/src/references.bib#popinet2009) used here
to test the compressible Navier--Stokes solver with VOF interface tracking and
surface tension. */

#define JACOBI 1

#include "grid/multigrid.h"
#include "two-phase-compressible.h"
#include "compressible-tension.h"

#define DIAMETER 0.8
#define MU sqrt(DIAMETER/LAPLACE)
#define TMAX (sq(DIAMETER)/MU)

/**
We will vary the number of levels of refinement (to study the
convergence) */

int LEVEL;
double LAPLACE;
double DC = 0.;
FILE * fp = NULL;

int main() {
    
  PI1 = 300.;
  gamma1 = 7.14;
  gamma2 = 1.4;

  TOLERANCE = 1e-4;
  f.sigma = 1;
  f.gradient = zero;

  /**
  Wne fix the Laplace number and look for stationary solutions
  (i.e. the volume fraction field is converged to within 1e-10) for
  varying spatial resolutions. */

  LAPLACE = 12000; DC = 1e-10;
  for (LEVEL = 3; LEVEL <= 6; LEVEL++) {
    N = 1 << LEVEL;
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

event logfile (i += 100; t <= TMAX)
{
  printf("%g \n", t/TMAX);
}

event error (t = end) {
  
  /**
  At the end of the simulation, we compute the equivalent radius of
  the droplet. */

  double vol = statsf(f).sum;
  double radius = sqrt(4.*vol/pi);

  /**
  And compute the maximum error on the curvature *ekmax*, the norm of
  the velocity *un* and the shape error *ec*. */


  double ekmax = 0., eshape = 0., perim = 0., eshapemax = 0;
  scalar un[], ec[], kappa[];
  curvature (f, kappa);
  foreach() {
    un[] = norm(q);

    /** 
    Integration of the error on the interface position with respect to a perfect circle. */
    face vector s;
    s.x.i = -1;
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = facet_normal (point, f, s);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2) {
        double length = sqrt(sq(segment[0].x*Delta - segment[1].x*Delta) 
                           + sq(segment[0].y*Delta - segment[1].y*Delta));
        double error1 = fabs(sqrt(sq(x + segment[0].x*Delta) + sq(y + segment[0].y*Delta)) - radius);
        double error2 = fabs(sqrt(sq(x + segment[1].x*Delta) + sq(y + segment[1].y*Delta)) - radius);
        
        perim  += length;
        eshape += (error1 + error2)/2.*length;
        eshapemax = max(eshapemax, max(error1,error2));
      }
    }

    /**
    We obtain also the error in curvature */
    if (kappa[] != nodata) {
      double ek = fabs (kappa[] - (/*AXI*/ + 1.)/radius);
      if (ek > ekmax)
        ekmax = ek;
    }

  }
  
  /**
  We output these on standard error (i.e. the *log* file). */

  fprintf (stderr, "%d %g %g %g %g %g %g\n", 
     LEVEL, LAPLACE, 
     normf(un).max*sqrt(DIAMETER), 
     perim, eshape/perim, eshapemax,
     ekmax);
}

/**
## Results

~~~gnuplot Convergence of the error on the equilibrium shape of the droplet with resolution. The diameter is given in number of grid points.
set xlabel 'D'
set ylabel 'Shape error'
set logscale xy
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):5 w lp t "RMS" ps 2, \
            '< sort -n -k1,2 /home/popinet/basilisk/src/test/spurious/log' u (0.8*2**$1):5 w lp t "RMS (incomp)" ps 2, \
            '< sort -n -k1,2 log' u (0.8*2**$1):6 w lp t "Max" ps 2, \
            '< sort -n -k1,2 /home/popinet/basilisk/src/test/spurious/log' u (0.8*2**$1):6 w lp t "Max (incomp)" ps 2, \
             0.2/(x*x) t "Second order"
~~~

~~~gnuplot Convergence of the relative error on the equilibrium curvature value with resolution. The diameter is given in number of grid points.
set ylabel 'Relative curvature error'
set log xy
plot [5:120]'< sort -n -k1,2 log' u (0.8*2**$1):($7/2.5) w lp t "Max" ps 2, \
             0.6/(x*x) t "Second order", \
             '< sort -n -k1,2 /home/popinet/basilisk/src/test/spurious/log' u (0.8*2**$1):($7/2.5) w lp t "Ref (incompressible)" ps 2
~~~

## See also

* [Effect of the Laplace number](http://basilisk.fr/sandbox/fuster/Allmach3.0/spurioustest.c)
* [Same test with the incompressible version](http://basilisk.fr/src/test/spurious.c)
* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/
