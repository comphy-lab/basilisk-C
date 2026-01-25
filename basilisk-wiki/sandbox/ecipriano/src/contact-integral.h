/**
# Contact angle from a distance function

The idea is to mimic the implementation of [contact.h](/src/contact.h) but in
simulations employing the [two-phase-clsvof.h](/src/two-phase-clsvof.h) with [integral.h](/src/integral.h) surface tension method.

The contact angle seems to be correct when looking at the distance function
isoline, but the VOF field creates strange artifacts in the vicinity of the
contact line. 

## Boundary condition

Using the following boundary condition for the distance field guarantees the
correct contact angle.
*/

#define contact_angle(theta)    \
  (val(_s) == nodata ? nodata : val(_s) + Delta*cos(theta))

/**
## Interface normals

We overwrite the interface normals calculation for VOF in order to account for
the modified interface position in the ghost cells. */

coord interface_normal (Point point, scalar c);

#undef interface_normal
#define interface_normal(point, c) interface_normal (point, c)

extern scalar d;

coord interface_normal (Point point, scalar c)
{
  double dx = (d[1] - d[-1])/2.;
  double dy = (d[0,1] - d[0,-1])/2.;
  double dn = sqrt(sq(dx) + sq(dy)) + 1e-30;
  coord n;
  n.x = dx / dn;
  n.y = dy / dn;
  return n;
}

