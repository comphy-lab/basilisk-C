/**
# Contact angles

This file is used to impose contact angles on boundaries for
interfaces described using a [VOF](vof.h) tracer and [height
functions](heights.h).

We first overload the default function used to compute the normal,
defined in [fractions.h](). */

#include "fractions.h"
#include "curvature.h"

/**
We will compute the normal using height-functions instead. If this is
not possible (typically at low resolutions) we revert back to
the Mixed-Youngs-Centered approximation. */

coord height_myc_normal (Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  return n;
}

macro coord interface_normal (Point point, scalar c) {
  return height_myc_normal (point, c);
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
angles. Follow the documentation in the section [Small/Large contact
angles](#smalllarge-contact-angles) to activate support for very small (or very
large) contact angle values. */

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
## Small/Large contact angles

The small contact angles implementation can be used by setting additional
boundary conditions, on the volume fraction and on the normal components
of the height functions:

~~~literatec
f[boundary] = contact_fraction (theta0*pi/180.);
h.n[boundary] = contact_normal (theta0*pi/180., f);
~~~

The implementation reflects the method proposed by [Afkhami \& Bussmann,
2009](#afkhami2008) and, as such, it is limited to 2D. */

/**
First, we define functions which allows to categorize the cells in the vicinity
of the contact line. In particular:

* *is_contact*: defines a cell along the boundary which contains the contact line
* *is_adjacent*: defines an interfacial cell next to the contact line cell along
  the boundary direction.
* *is_opposite*: defines a non-interfacial cell next to the contact line cell
  along the boundary direction.
*/

foreach_dimension()
static inline int is_contact_x (Point point, scalar c, double theta0) {
  if (!is_leaf (cell))
    return false;
  if (c[] == 0. || c[] == 1.)
    return false;
  else if (is_boundary (neighbor(1)) || is_boundary (neighbor(-1))) {
      for (int i = -1; i <= 1; i += 2)
        if (theta0 < pi/2. && c[0,i] == 0.)
          return true;
        else if (theta0 >= pi/2. && c[0,i] == 1.)
          return true;
      }
  return false;
}

foreach_dimension()
static inline int is_adjacent_x (Point point, scalar c, double theta0) {
  if (!is_leaf (cell))
    return false;
  if (c[] == 0. || c[] == 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (is_contact_x (neighborp(0,i), c, theta0))
        return true;
  }
  return false;
}

foreach_dimension()
static inline int is_opposite_x (Point point, scalar c, double theta0) {
  if (!is_leaf (cell))
    return false;
  if (c[] > 0. && c[] < 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (is_contact_x (neighborp(0,i), c, theta0)) {
        if (theta0 < pi/2. && c[] == 0.)
          return true;
        else if (theta0 >= pi/2. && c[] == 1.)
          return true;
      }
  }
  return false;
}

/**
We compute the normal pointing outward from the domain boundary. */

static inline
coord normal_boundary (Point point, Point neighbor) {
  return (coord){neighbor.i - point.i,
                 neighbor.j - point.j};
}

/**
The following function returns the normal of an interface touching the boundary.
`nb` is the normal to the boundary, `nf` is the VOF interface normal (without
considering the contact angle), while `angle` is the contact angle value. */

static inline
coord normal_contact (coord ns, coord nf, double angle)
{
  coord n;
  if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
    n.x = - ns.x*cos(angle) + ns.y*sin(angle);
    n.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else {
    n.x = - ns.x*cos(angle) - ns.y*sin(angle);
    n.y =   ns.x*sin(angle) - ns.y*cos(angle);
  }
  return n;
}

/**
### Volume fraction boundary condition

We adjust the value of the volume fraction in the ghost cells in order to match
the interface orientation promoted by the contact angle. */

foreach_dimension()
static double contact_fraction_x (double expr, Point point = point,
    Point neighbor = neighbor, scalar s = _s)
{
  if (is_contact_x (point, s, expr)) {
    coord m = mycs (point, s);
    coord mb = normal_boundary (point, neighbor);
    coord mcl = normal_contact (mb, m, expr);
    double alpha = plane_alpha (s[], mcl);

    coord a = {-0.5, -0.5}, b = {0.5, 0.5};
    foreach_dimension() {
      a.x += mb.x;
      b.x += mb.x;
    }
    return rectangle_fraction (mcl, alpha, a, b);
  }
  else if (s[] > 0.)
    return expr < pi/2. ? 1. : 0.;
  else
    return s[];
}

double contact_fraction (double expr, Point point = point,
    Point neighbor = neighbor, scalar s = _s)
{
  if (neighbor.i != point.i)
    return contact_fraction_x (expr, point, neighbor, s);
  if (neighbor.j != point.j)
    return contact_fraction_y (expr, point, neighbor, s);
  assert (false);
  return 0;
}

/**
### Normal component of the height functions

This boundary conditions works in an unusual manner: instead of setting the
ghost value of the height function, it directly imposes the right value in the
contact line and the opposite cells. This is required because the curvature
computed using the normal component of the height function along the boundary
does not call the ghost cell values, and it is a peculiarity of the small/large
angle cases. */

foreach_dimension()
static inline
double height_contact_x (Point point, scalar c, double expr, coord mb) {
  coord m = mycs (point, c);
  coord mcl = normal_contact (mb, m, expr);
  double alpha = plane_alpha (c[], mcl);
  double hc = nodata;
  if (fabs (mcl.x) > 0. && fabs (alpha/mcl.x) <= 5.5)
    hc = alpha/mcl.x + (mcl.x < 0.)*HSHIFT;
  return hc;
}

foreach_dimension()
static double contact_normal_x (double expr, scalar c,
    Point point = point, Point neighbor = neighbor, scalar s = _s)
{
  if (is_contact_x (point, c, expr)) {
    coord mb = normal_boundary (point, neighbor);
    double hc = height_contact_x (point, c, expr, mb);
    if (s[] == nodata)
      s[] = hc;
    return s[];
  }
  if (is_opposite_x (point, c, expr)) {
    int i = is_contact_x (neighborp(0,1), c, expr) ? 1 :
            is_contact_x (neighborp(0,-1), c, expr) ? -1 : 0;
    if (i == 0)
      return nodata;
    else {
      coord mb = normal_boundary (point, neighbor);
      double hc = height_contact_x (neighborp(0,i), c, expr, mb);
      if (hc != nodata) {
        s[] = (point.level < neighborp(0,i).level) ?
          hc + (orientation (hc) ? 1. : -1.)*1.5*tan (expr) :
          hc + (orientation (hc) ? 1. : -1.)*tan (expr);
        return s[];
      }
    }
  }
  return nodata;
}

double contact_normal (double expr, scalar c,
    Point point = point, Point neighbor = neighbor, scalar s = _s)
{
  if (neighbor.i != point.i)
    return contact_normal_x (expr, c, point, neighbor, s);
  if (neighbor.j != point.j)
    return contact_normal_y (expr, c, point, neighbor, s);
  assert (false);
  return 0;
}

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

@article{afkhami2008,
  title={Height functions for applying contact angles to 2D VOF simulations},
  author={Afkhami, Shahriar and Bussmann, Markus},
  journal={International journal for numerical methods in fluids},
  volume={57},
  number={4},
  pages={453--472},
  year={2008},
  publisher={Wiley Online Library},
  url={https://web.njit.edu/~shahriar/Publication/IJNMF1.pdf}
}
~~~
*/
