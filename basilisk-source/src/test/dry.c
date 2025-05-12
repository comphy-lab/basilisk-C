/**
# Drying of a lake

We test the robustness of the explicit and implicit Saint-Venant
solvers for an emptying lake in one dimension. The topography dries
and the lake-at-rest condition must be reached at the end of the
run. */

#include "grid/multigrid1D.h"
#if EXPLICIT
# include "saint-venant.h"
#else
# include "saint-venant-implicit.h"
#endif

int main()
{
  origin (-0.5);
  init_grid (256);
  DT = 1e-2;
  run();
}

#if EXPLICIT
h[right] = dry;
eta[right] = zb[];
u.n[right] = u.x[];
#else
q.n[right] = q.x[];
p[right] = dirichlet(0);
#endif

event init (i = 0)
{
  double a = 0.25, k = pi/0.1;
  foreach() {
    zb[] = a*(cos(k*x) + 1.)*(fabs(x) < 0.1);
    h[] = 0.8 - zb[];
  }
}

event logfile (t = {0.5, 0.75, 1, 3, 50}) {
  printf ("%g %g %d\n", t, dt, i);
  foreach() {
    fprintf (stderr, "%g %g %.6f %g\n", x, h[],
#if EXPLICIT
	     u.x[],
#else
	     q.x[],
#endif
	     zb[]);
    assert (h[] > 0.);
  }
  fprintf (stderr, "\n");
}

/**
~~~gnuplot Evolution of the free surface (red) and topography (green)
set term @SVG enhanced size 640,200
unset key
plot 'log' u 1:($2+$4) w l, 'log' u 1:4 w l
~~~
*/
