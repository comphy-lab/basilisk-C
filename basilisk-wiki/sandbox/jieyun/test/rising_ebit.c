/**
# Rising bubble

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can solve either the
axisymmetric or planar version. We can used standard or "reduced"
gravity. */

#define LEVEL 8
#define ADAPT 1

face vector av[];

#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif

#include "navier-stokes/centered.h"
#include "two-phase-ebit.h"
#include "tension.h"

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);

/**
We make sure there is no flow through the top and bottom boundary,
otherwise the compatibility condition for the Poisson equation can be
violated. */

uf.n[bottom] = 0.;
uf.n[top] = 0.;

FILE * fp_dis;

int main() {
  TOLERANCE = 1e-4 [*];
  DT = 1.e-2 [0, 1];
  rho1 = 1000. [0];
  mu1 = 10.;
  double sigma_in;

  #ifdef CASE2
  rho2 = 1.;
  mu2 = 0.1;
  sigma_in = 1.96;
  TOLERANCE = 1e-3 [*];
  DT = 1.e-3;
  #else
  rho2 = 100.;
  mu2 = 1.;
  sigma_in = 24.5;
  #endif

  f.sigma = sigma_in;

  size (2. [1]);
  init_grid (1 << LEVEL);
  a = av;

  run();
  fclose (fp_dis);
}


event init (i = 0) {
  /**
  The domain is a rectangle. We only simulate half the bubble. */
  
  mask (y > 0.5 ? top : none);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x - 0.5) + sq(y) - sq(0.25);
  
  init_markers (phi);

  char name[80];
  sprintf (name, "rising_dis.dat");
  fp_dis = fopen (name, "w");
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  foreach_face(x)
    av.x[] -= 0.98;
  boundary ((scalar *){av});
}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({mask_intf, u}, (double[]){0.02, 1e-3, 1e-3}, maxlevel = LEVEL, minlevel = LEVEL - 3);
}
#endif

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
  }

  fprintf (fp_dis, "%g %g %g %g %g %g %g %g \n", 
	  t, sb, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
  fflush (fp_dis);
}


/**
At $t=3$ we output the shape of the bubble. */

event interface (t = 3.) {
  char out[100];
  sprintf (out, "rising_intf_ebit.dat");
  output_facets_ebit (out);
}

/**
## Results

The final shape of the bubble is compared to that obtained with the VOF method in the Basilisk
[Rising bubble](http://basilisk.fr/src/test/rising.c) test.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set term push
set term @SVG size 640,320
set size ratio -1
set grid
plot [][0:0.4]'http://basilisk.fr/src/test/rising/log' u 1:2 w l t 'Basilisk', \
              'rising_intf_ebit.dat' u 1:2 w l t 'EBIT'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'http://basilisk.fr/src/test/rising/out' u 1:5 w l t 'Basilisk', \
              'rising_dis.dat' u 1:5 w l t 'EBIT'
~~~

*/
