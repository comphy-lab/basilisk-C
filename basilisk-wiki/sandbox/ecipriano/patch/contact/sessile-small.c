/**
# Sessile drop with small contact angles

This test case is similar to [sessile.c](/src/test/sessile.c), and it demonstrates
the accuracy of the additional boundary conditions for very small (or very large)
contact angle values.

~~~gnuplot Equilibrium shapes
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-2:2]
set yrange [0:]
plot 'out' u 2:1 w l, '' u (-$2):1 w l lt 1, 0 lt -1
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact-small.h"
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};

/**
To set the contact angle, we allocate a [height-function field](/src/heights.h)
and set the contact angle boundary condition on its tangential component.
*/

vector h[];
double theta0 = 30;
h.t[left] = contact_angle (theta0*pi/180.);

/**
We activate the support for small angles using the following additional boundary
conditions. */

f[left] = contact_fraction (theta0*pi/180.);
h.n[left] = contact_normal (theta0*pi/180., f);

int LEVEL;

int main()
{
  size (3 [1]);

  /**
  We use a constant viscosity. */

  mu[] = {.1,.1};

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;


  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */

  f.sigma = 1.;

  for (LEVEL = 4; LEVEL <= 6; LEVEL++) {
    for (theta0 = 10.; theta0 <= 30.; theta0 += 10.) {
      init_grid (1 << LEVEL);
      run();
    }
  }
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, {1e-3,1e-2,1e-2}, LEVEL, 3);
}
#endif

/**
At equilibrium (t = 20 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

event end (t = 30)
{
  if (LEVEL == 6)
    output_facets (f, stdout);

  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "level %d %g %.5g %.3g\n", LEVEL, theta0, R/sqrt(V/pi), s.stddev);
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set xr[10:30]
set xtics 10,10,30
set grid
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  '<grep "level 4" log' u 3:4 pt 7 t 'LEVEL 4', \
  '<grep "level 5" log' u 3:4 pt 7 t 'LEVEL 5', \
  '<grep "level 6" log' u 3:4 pt 7 t 'LEVEL 6'
~~~
*/
