/**
# Pinning a Liquid Droplet

We want to use the method proposed in [sandbox/lopez/contact.c](/sandbox/lopez/contact.c)
to pin the liquid droplet at a specific point of the domain.

This configuration is important, from the experimental
point of view, to suspend droplets on a solid fiber
in normal gravity conditions, in order to study the
evaporation and combustion characteristics of that
droplet under the influence of gravity.

We select a point where the droplet must be pinned,
and we let the gravity field deform the liquid droplet
until reaching a steady position.

~~~gnuplot Droplet shape at different gravity values
set term push
set term @SVG size 640,180
set size ratio -1
#unset key
set xrange[-1.5:1.5]
set yrange[0:0.7]
unset xtics
unset ytics
unset border
plot 0 lt -1 notitle, \
     "out-0.0" w l t "G.x = -0.0", \
     "out-0.5" w l t "G.x = -0.5", \
     "out-1.0" w l t "G.x = -1.0", \
     "out-1.5" w l t "G.x = -1.5", \
     "out-2.0" w l t "G.x = -2.0", \
     "out-2.5" w l t "G.x = -2.5"
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

vector h[];

/**
Known the position of the interface $a$ we compute the corresponding
heigth function at the neighbouring cell.*/

foreach_dimension()
static double line_x (Point point, scalar h, double a)
{
  return (h[] == nodata ? nodata : -h[] + 2.*(a-x)/Delta);
}

#define contact_line(theta) contact_line_ (point, neighbor, _s, theta)

double contact_line_ (Point point, Point neighbor, scalar h, double a)
{
  if (neighbor.i != point.i)
    return line_x (point, h, a);
  if (neighbor.j != point.j)
    return line_y (point, h, a);
  assert (false); // not reached
  return 0.;
}

double ap;
h.t[bottom] = x > 0 ? contact_line (ap) : neumann (0);

int main()
{
  size (2);
  origin (-1);

  /**
  We set a small density and viscosity ratio. */

  mu1 = 1.;
  mu2 = 0.1;
  rho2 = 0.1;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set the surface tension coefficient and run for two position of
  the pinning point. */

  f.sigma = 1;
  ap = 0.75;

  /**
  We run the simulation at different values of gravity .*/

  for (G.x = 0.; G.x >= -2.5; G.x -= 0.5)
    run();
}

/**
We initialize half liquid droplet on the bottom boundary. */

event init (t = 0)
{
  fraction (f, - (sq(x - 0.25) + sq(y) - sq(0.5)));
}

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape. */

event end (t = 10)
{
  char name[80];
  sprintf (name, "out-%.1f", fabs(G.x));

  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

