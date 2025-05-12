/**
# Anchorage of the contact line

A sessile droplet is a drop of liquid at rest on a solid surface. In
the absence of gravity, the shape of the droplet is controlled only by
the surface tension and the particular liquid-solid interaction.
Normally, the  [contact
angle](/src/contact.h) or the position of the contact line can be
imposed. Here we propose a method to force the position of the contact
line.

To test the proposed procedure, a droplet that is initialized as a
half disk is forced to anchor at different positions. The droplet
oscillates and finally relaxes to its equilibrium position.

Note that if the desired anchoring position is too far away from the
initial position of the interface, the interface will not be attracted
to the desired position.

~~~gnuplot Equilibrium shapes for contact line at $a$=0.4 and $a$=0.6.
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.6:1.6]
set yrange [0:]
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};
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
h.t[bottom] = contact_line (ap);

int main()
{
  size (2);

  /**
  We use a constant viscosity. */
  
  const face vector muc[] = {.1,.1};
  mu = muc;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set the surface tension coefficient and run for two position of
  the contact line */
  
  f.sigma = 1.;

  ap=0.4;
  run();
  
  ap=0.6;
  run();
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape. */

event end (t = 10)
{
  output_facets (f, stdout);
}

