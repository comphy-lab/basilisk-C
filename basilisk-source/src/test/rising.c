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
gravity. We also test levelset interface tracking and a momentum
formulation. */

#include "grid/multigrid.h"
#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif
#if MOMENTUM
# include "momentum.h"
#else
#include "navier-stokes/centered.h"
#if CASE2
# define FILTERED 1
#endif
#if CLSVOF
# include "two-phase-clsvof.h"
#elif LEVELSET
# include "two-phase-levelset.h"
#else
# include "two-phase.h"
#endif
#endif
#if LEVELSET
# include "integral.h"
#else
# include "tension.h"
#endif
#if REDUCED
# include "reduced.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

#if MOMENTUM
q.t[right] = dirichlet(0);
q.t[left]  = dirichlet(0);
#else
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
#endif

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  dimensions (nx = 4);
  size (2 [1]);
  DT = 1. [0,1];
  init_grid (1 << LEVEL);
  
  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */

  rho1 = 1000.[0], mu1 = 10.;  // works also with rho1 = 1000.[-3,0,1]
#if CASE2
  rho2 = 1., mu2 = 0.1;
#else
  rho2 = 100., mu2 = 1.;
#endif

#if LEVELSET
  #if CASE2
  const scalar sigma[] = 1.96;
  #else
  const scalar sigma[] = 24.5;
  #endif
  d.sigmaf = sigma;
#else // !LEVELSET
  #if CASE2
  f.sigma = 1.96;
  #else
  f.sigma = 24.5;
  #endif
#endif // !LEVELSET
  
  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4 [*];
#if REDUCED
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

event init (t = 0) {

  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */

#if LEVELSET
  foreach()
    d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
#else
  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
#endif
}

/**
We add the acceleration of gravity. */

#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}
#endif

/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
#if MOMENTUM
    vb += q.x[]*dv/rho(f[]);
#else
    vb += u.x[]*dv;
#endif
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    printf ("t sb -1 xb vb dt perf.t perf.speed\n");
    sb0 = sb;
  }
  printf ("%g %g %g %g %g %g %g %g ", 
	  t, (sb - sb0)/sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
#if !MOMENTUM
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
#endif
  putchar ('\n');
  fflush (stdout);
}

/**
At $t=3$ we output the shape of the bubble. */

event interface (t = 3.) {
  output_facets (f, stderr);
}

#if _GPU && SHOW
event display (i++) {
  output_ppm (f, fp = NULL, fps = 30, spread = -1, n = 1024);
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, fp = NULL, fps = 30, spread = 6, n = 1024);
}
#endif

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
}
#endif

/**
## Results

The final shape of the bubble is compared to that obtained with the
MooNMD Lagrangian solver (see [the FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html))
at the highest resolution. We also display the shape of the
axisymmetric version of the test. The axisymmetric bubble moves much
faster.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set term push
set term @SVG size 640,320
set size ratio -1
set grid
plot [][0:0.4]'../c1g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              'log' u 1:2 w l t 'Basilisk', \
              '../rising-levelset/log' u 1:2 w l t 'Basilisk (levelset)', \
              '../rising-clsvof/log' u 1:2 w l t 'Basilisk (CLSVOF)', \
              '../rising-axi/log' u 1:2 w l t 'Basilisk (axisymmetric)', \
              '../rising-axi-clsvof/log' u 1:2 w l t 'Basilisk (axi + CLSVOF)', \
              '../rising-axi-momentum/log' u 1:2 w l t 'Basilisk (axi + momentum)'
~~~

For test case 2, the mesh in Basilisk is too coarse to accurately
resolve the skirt.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 2.
set key bottom left
plot [][0:0.4]'../c2g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              '../rising2/log' u 1:2 w l t 'Basilisk', \
              '../rising2-levelset/log' u 1:2 w l t 'Basilisk (levelset)', \
              '../rising2-clsvof/log' u 1:2 w l t 'Basilisk (CLSVOF)'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c1g3l4.txt' u 1:5 w l t 'MooNMD', \
              'out' u 1:5 w l t 'Basilisk', \
              '../rising-levelset/out' u 1:5 w l t 'Basilisk (levelset)', \
              '../rising-clsvof/out' u 1:5 w l t 'Basilisk (CLSVOF)',     \
              '../rising-axi/out' u 1:5 w l t 'Basilisk (axisymmetric)',  \
              '../rising-axi-clsvof/out' u 1:5 w l t 'Basilisk (axi + CLSVOF)',  \
              '../rising-axi-momentum/out' u 1:5 w l t 'Basilisk (axi + momentum)'
~~~

~~~gnuplot Relative volume difference as a function of time for test case 1.
reset
set grid
set xlabel 'Time'
set ylabel '(vb - vb_0)/vb_0'
set key bottom left
plot [0:3]'out' u 1:2 w l t 'Basilisk', \
          '../rising-levelset/out' u 1:2 w l t 'Basilisk (levelset)',   \
	  '../rising-clsvof/out' u 1:2 w l t 'Basilisk (CLSVOF)',	\
	  '../rising-axi/out' u 1:2 w l t 'Basilisk (axisymmetric)',	\
	  '../rising-axi-clsvof/out' u 1:2 w l t 'Basilisk (axi + CLSVOF)',	\
	  '../rising-axi-momentum/out' u 1:2 w l t 'Basilisk (axi + momentum)'
~~~

~~~gnuplot Rise velocity as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c2g3l4.txt' u 1:5 w l t 'MooNMD', \
              '../rising2/out' u 1:5 w l t 'Basilisk', \
              '../rising2-levelset/out' u 1:5 w l t 'Basilisk (levelset)', \
              '../rising2-clsvof/out' u 1:5 w l t 'Basilisk (CLSVOF)'
~~~

~~~gnuplot Relative volume difference as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set ylabel '(vb - vb_0)/vb_0'
set key top left
plot [0:3]'../rising2/out' u 1:2 w l t 'Basilisk',		       \
          '../rising2-levelset/out' u 1:2 w l t 'Basilisk (levelset)', \
	  '../rising2-clsvof/out' u 1:2 w l t 'Basilisk (CLSVOF)'
~~~
*/
