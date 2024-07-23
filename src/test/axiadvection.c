/**
# Axisymmetric mass conservation

A standard and a VOF tracer are advected by an axisymmetric flow. The
initial interface is a torus which is then advected by the flow
illustrated in the figure below. As the torus is flattened against the
right-hand-side wall, its cross-sectional surface area decreases but
the volume should remain constant. 

~~~gnuplot Evolution of the VOF interface and velocity field
set size ratio -1
set xlabel 'z'
set ylabel 'r'
plot [-0.5:0.5][0:1]'out' w l t '', 'velo' u 1:2:($3/17.):($4/17.) w vect t ''
~~~
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tracer.h"

scalar f[], f1[];
scalar * interfaces = {f}, * tracers = {f1};

int main()
{
  X0 = -0.5;
  N = 64;
  TOLERANCE = 1e-12 [*];
  f1.gradient = minmod2;
  run();
}

u.n[left] = dirichlet(1);
u.t[left] = dirichlet(0);
p[left]   = neumann(0);

u.n[top] = neumann(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);

#define ellipse(xc, yc, a, b) (sq((x - xc)/(a)) + sq((y - yc)/(b)) - 1.)

event init (i = 0) {
  foreach()
    u.x[] = 1.;
  fraction (f, - ellipse (0, 0.3, 0.1, 0.1));
  fraction (f1, - ellipse (0, 0.3, 0.1, 0.1));
}

event logfile (i++; t <= 0.8) {
  static double sfmin = HUGE, sfmax = -HUGE;
  static double sfmin1 = HUGE, sfmax1 = -HUGE;
  double s = statsf(f).sum, s1 = statsf(f1).sum;
  if (s < sfmin) sfmin = s;
  if (s > sfmax) sfmax = s;
  if (s1 < sfmin1) sfmin1 = s1;
  if (s1 > sfmax1) sfmax1 = s1;
  double e = 2.*(sfmax - sfmin)/(sfmax + sfmin);
  double e1 = 2.*(sfmax1 - sfmin1)/(sfmax1 + sfmin1);
  fprintf (stderr, "%g %.12f %.12f %.10f %.10f\n", t, s, s1, e, e1);
  fflush (stderr);
  assert (e < 4e-8);
  assert (e1 < 5e-5);
}

event output (t += 0.2; t <= 1.2)
  output_facets (f);

event velo (t = end)
  output_field ((scalar *){u}, fopen ("velo", "w"), n = 16, linear = true);

#if TREE

#if 0
event gfsview (i++) {
  static FILE * fp = popen ("gfsview2D -s test.gfv", "w");
  output_gfs (fp);
}
#endif

event adapt (i++) {
  double sb = statsf(f).sum;
  double sb1 = statsf(f1).sum;
  adapt_wavelet ({f1}, (double[]){5e-3}, maxlevel = 6, minlevel = 4);
  double sa = statsf(f).sum;
  double sa1 = statsf(f1).sum;
  // the mass of VOF tracers is not conserved exactly
  assert (fabs(sa - sb) < 2e-6);
  // the mass of diffusive tracers must be conserved to within round-off
  assert (fabs(sa1 - sb1) < 1e-12);
}
#endif

/**
## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/axiadvection.html)
*/
