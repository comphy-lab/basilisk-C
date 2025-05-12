/**
# Rising bubble conserving

Modified version of the original code with conserving and without the famous mask function..
Comparison with and without adapt case 1 and 2 are done with a patch concerning vof-tracer not yet in /src.

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can solve either the
axisymmetric or planar version. We can used standard or "reduced"
gravity. */

#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif
#include "navier-stokes/centered.h"
#include "two-phase.h"
#if CONSERVING
#include "navier-stokes/conserving.h"
#endif
#include "tension.h"
#if REDUCED
# include "reduced.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif

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

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  size (2);
  init_grid (1 << LEVEL);
  
  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */

  rho1 = 1000., mu1 = 10.;
#if CASE2
  rho2 = 1., mu2 = 0.1, f.sigma = 1.96;
#else
  rho2 = 100., mu2 = 1., f.sigma = 24.5;
#endif

  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4;
#if REDUCED
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

event init (t = 0) {

  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */

  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
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
    xb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
  }
  printf ("%g %g %g %g %g %g %g %g ", 
	  t, sb, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}

/**
At $t=3$ we output the shape of the bubble. */

event interface (t = 3.) {
  output_facets (f, stderr);
}

/**
If bview is installed on the system, we can also visualise the
simulation as it proceeds. */

#if 0
event bview (i += 10) {
  char name[80];
  scalar omega[];
  vorticity (u, omega);
  sprintf (name, "snapshot-%g", t);
  dump (name);
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
at the highest resolution.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set term push
set term @SVG size 640,320
set size ratio -1
set grid
set key out right
plot [][0:0.5]'logn' u 1:2 w l t 'Basilisk',\
'log' u 1:2 w l t 'Basilisk with conserving',\
'logpatch' u 1:2 w l t 'Basilisk with conserving and patch',\
'logpatchadapt' u 1:2 w l t 'same with adapt mesh'
~~~

For test case 2, the mesh in Basilisk is too coarse to accurately
resolve the skirt. 
The test is not passed with conserving and without the patch.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 2.
plot [][0:0.5]'log2n' u 1:2 w l t 'Basilisk',\
'log2' u 1:2 w l t 'Basilisk with conserving',\
'log2patch' u 1:2 w l t 'Basilisk with conserving and patch',\
'log2patchadapt' u 1:2 w l t 'same with adapt mesh'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'outn' u 1:5 w l t 'Basilisk',\
'out' u 1:5 w l t 'Basilisk with conserving',\
'outpatch' u 1:5 w l t 'Basilisk with conserving and patch',\
'outpatchadapt' u 1:5 w l t 'same with adapt mesh'
~~~

The test is not passed with conserving and without the patch.

~~~gnuplot Rise velocity as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'out2n' u 1:5 w l t 'Basilisk',\
'out2' u 1:5 w l t 'Basilisk with conserving',\
'out2patch' u 1:5 w l t 'Basilisk with conserving and patch',\
'out2patchadapt' u 1:5 w l t 'same with adapt mesh'
~~~

Here we output the shape and min an max vertical velocity value at the last time step of the simulation.

![VOF representation with maximum and minimum value for vertical velocity colored, case 2 without conserving](rising-conserving/rising2.png)

![VOF representation with maximum and minimum value for vertical velocity colored, case 2 with conserving and without patch](rising-conserving/rising2-conserving.png)

![VOF representation with maximum and minimum value for vertical velocity colored, case 2 with conserving and with patch](rising-conserving/rising2p-conserving.png)

*/
