/**
# Pin a Liquid Droplet at a Specific Point

This module defines the boundary conditions required to
pin the droplet a specific point of the domain. The method
implemented here was extended from [sandbox/lopez/contact.c](/sandbox/lopez/contact.c).
*/

#include "contact.h"

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

/**
We create the local height-function vector and set the
boundary conditions. The values of ap and ac must be
tuned by the user:

* *ap*: x coordinate of the pinning point
* *ac*: x coordinate of the droplet centroid
*/

struct Pinning {
  double ap, ac;
};

struct Pinning pinning = {
  .ap = 0.,
  .ac = 0.,
};

vector h[];
h.t[bottom] = x > pinning.ac ? contact_line (pinning.ap) : neumann (0.);

/**
We must associate the height function field with the VOF tracer, so
that it is used by the relevant functions (curvature calculation in
particular). */

extern scalar f;

event defaults (i = 0) {
  f.height = h;
#ifdef AXI
  if (Y0 == 0.) {
    fprintf (ferr,
        "WARNING: Setting the contact line on the axis of symmetry, shift the origin along y\n");
    fflush (ferr);
  }
#endif
}

