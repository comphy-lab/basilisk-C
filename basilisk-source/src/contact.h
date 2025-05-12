/**
# Contact angles

This file is used to impose contact angles on boundaries for
interfaces described using a [VOF](vof.h) tracer and [height
functions](heights.h).

We first overload the default function used to compute the normal,
defined in [fractions.h](). */

coord interface_normal (Point point, scalar c);

#undef interface_normal
#define interface_normal(point, c) interface_normal (point, c)

#include "fractions.h"
#include "curvature.h"

/**
We will compute the normal using height-functions instead. If this is
not possible (typically at low resolutions) we revert back to
the Mixed-Youngs-Centered approximation. */

coord interface_normal (Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  return n;
}

/**
The height functions are stored in the vector field associated with
each VOF tracer. They need to be updated every time the VOF field
changes. For the [centered Navier-Stokes
solver](navier-stokes/centered.h), this means after initialisation and
after VOF advection. 

Note that strictly speaking this should be done for each
[sweep](vof.h#sweep_x) of the direction-split VOF advection, which we
do not do here i.e. we use the normal at the beginning of the timestep
and assume it is constant during each sweep. This seems to work
fine. */

extern scalar * interfaces;

event init (i = 0) {
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights (c, c.height);
}

event vof (i++) {
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights (c, c.height);
}

/**
The macro below can be used to impose a contact angle on a boundary by
setting the corresponding tangential component of the height
function. 

Note that the equivalent function for the normal component of the
height function is not defined yet. This limits the range of
accessible contact angles, since values of the normal component of the
height function will be required to compute curvature at shallow
angles. */

#if dimension == 2

#define contact_angle(theta)					\
  (val(_s) == nodata ? nodata : val(_s) +			\
   (orientation(val(_s)) ? -1. : 1.)/tan(theta))

/**
## Three-dimensional implementation

While the 2D implementation is trivial, in 3D one must take into
account the projection onto the boundary of the normal to the
interface (see [Afkhami \& Bussmann, 2009](#afkhami2009) for
details). This leads to the code below, where the only complication
comes from taking into account the relative orientations of the
boundary and height-function components. 

From a user point-of-view, using the *contact_angle()* macro is as
simple as in 2D. */

#else // dimension == 3

#define contact_angle(theta) contact_angle_ (point, neighbor, _s, theta)

foreach_dimension()
static double contact_z (Point point, scalar h, double theta)
{
  if (h.i == h.v.z.i) {
    fprintf (stderr,
	     "contact_angle() cannot be used for '%s' which is the normal\n"
	     "  component of the height vector\n",
	     h.name);
    exit (1);
  }

  if (h[] == nodata)
    return nodata;
  foreach_dimension(2)
    if (h.i == h.v.x.i)
      foreach_dimension(2) {
	coord n = normal2_x (point, h.v);
	if (n.x != nodata && n.y != nodata)
	  return h[] + 1./(tan(theta)*n.x/sqrt(sq(n.x) + sq(n.y)));
      }
  return h[]; // 90 degree contact angle if the normal is not defined
}

double contact_angle_ (Point point, Point neighbor, scalar h, double theta)
{
  if (neighbor.i != point.i)
    return contact_x (point, h, theta);
  if (neighbor.j != point.j)
    return contact_y (point, h, theta);
  if (neighbor.k != point.k)
    return contact_z (point, h, theta);
  assert (false); // not reached
  return 0.;
}

#endif // dimension == 3

/**
## References

~~~bib
@article{afkhami2009,
  title={Height functions for applying contact angles to 3D VOF simulations},
  author={Afkhami, S and Bussmann, M},
  journal={International Journal for Numerical Methods in Fluids},
  volume={61},
  number={8},
  pages={827--847},
  year={2009},
  publisher={Wiley Online Library},
  url={https://web.njit.edu/~shahriar/Publication/IJNMF2.pdf}
}
~~~
*/
