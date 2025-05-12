/**
# Evolution by mean curvature flow

Solve
$$
\mathbf{n}\cdot\partial_\tau\mathbf{X} = -\kappa
$$
See e.g [Deckelnick et al,
2005](https://homepages.warwick.ac.uk/~masgat/PAPERS/DecDziEll05_ActaNumerica.pdf).

~~~gnuplot Evolution of the interface
set size ratio -1
plot 'out' w l t ''
~~~
*/

#include "grid/multigrid.h"
#include "advection.h"
#include "vof.h"
#include "curvature.h"

scalar f[];
scalar * tracers = NULL, * interfaces = {f};

int main()
{
  L0 = 1 [0];
  origin (-0.5, -0.5);
  N = 64;
  // this is the diffusive timestep limit (with diffusion coeff unity)
  DT = sq(L0/N)/2.;
  run();
}

double circle (double x, double y)
{
  double theta = atan2(y,x);
  double r = sqrt(sq(x) + sq(y));
  return r - 0.4*(1. + 0.2*cos(8.*theta));
}

event init (i = 0)
{
  fraction (f, circle(x, y));
}

coord normal (Point point, scalar c)
{
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

scalar kappa[];

event velocity (i++) {
  curvature (f, kappa);
  foreach_face() {
    coord n;
    foreach_dimension()
      n.x = 0.;
    double k = 0.;
    if (f[] > 0. && f[] < 1.) {
      n = normal (point, f);
      k = kappa[];
      if (f[-1] > 0. && f[-1] < 1.) {
	coord np = normal (neighborp(-1), f);
	foreach_dimension()
	  n.x = (n.x + np.x)/2.;
	k = (k + kappa[-1])/2.;
      }
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      k = kappa[-1];
      n = normal (neighborp(-1), f);
    }
    u.x[] = - (k /*+ 2.5*/)*n.x;
  }
#if 0  
  FILE * fp = fopen("velo", "w");
  foreach()
    fprintf (fp, "%g %g %g %g\n", x, y,
	     (uf.x[] + uf.x[1])/2.,
	     (uf.y[] + uf.y[0,1])/2.);
  fclose (fp);
#endif
}

event logfile (t += 0.001; t <= sq(0.1)) {
  //  stats s = statsf (kappa);
  //  fprintf (stderr, "# %f %d %.12f %g %g\n", t, i, s.sum, s.min, s.max);
  output_facets (f, stdout);
}

#if 0
event gfsview (i++) {
  static FILE * fp = popen ("gfsview2D -s mean.gfv", "w");
  output_gfs (fp, t = t);
}
#endif
