/**
# Contact angles

This file is used to impose contact angles on boundaries for
interfaces described using a [VOF](vof.h) tracer and [height
functions](heights.h).

We first overload the default function used to compute the normal,
defined in [fractions.h](). */

#include "utils_ebm.h"

#include "fractions.h"
#include "curvature.h"

/**
We will compute the normal using height-functions instead. If this is
not possible (typically at low resolutions) we revert back to
the Mixed-Youngs-Centered approximation. */

coord height_myc_normal (Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata) {
    n = mycs (point, c);
    fprintf (stdout, "t = %g, height_normal fails at (x = %g, y = %g)!\n", t, x, y);
  }
  return n;
}

macro coord interface_normal (Point point, scalar c) {
  return height_myc_normal (point, c);
}

extern scalar * interfaces;

event init (i = 0) {
  for (scalar c in interfaces)
    if (c.height.x.i) {
#if !EMBED
      heights (c, c.height);
#else
      heights_ebm (c, cs, c.height, c.cet);
#endif
    }
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
    fprintf (stdout,
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


// 
bool is_three_phase (double c, double cs)
{
  if (c > 0. && c < cs && cs < 1.)
    return true;
  else
    return false;
}

coord find_cl_pos (coord c_n, double c_alpha, coord cs_n, double cs_alpha)
{
  coord cross;
  assert (fabs(cs_n.y*c_n.x - cs_n.x*c_n.y) > 0);
  cross.x = (c_alpha*cs_n.y - cs_alpha*c_n.y)/(cs_n.y*c_n.x - cs_n.x*c_n.y);
  assert (fabs(c_n.y) > 0 || fabs(cs_n.y) > 0);
  if (fabs(c_n.y) > fabs(cs_n.y))
    cross.y = (c_alpha - c_n.x*cross.x)/c_n.y;
  else
    cross.y = (cs_alpha - cs_n.x*cross.x)/cs_n.y;

  return cross;
}

#if dimension == 2

static void check_cfg (const scalar c, const scalar cs, face vector fs, scalar cfg)
{
  foreach()
    if (cs[] > 0. && cs[] < 1.) {
      assert (cfg[] != nodata);
      if (cfg[] == 1. && !is_three_phase (c[], cs[])) {
        fprintf (ferr, "CL locates in non-interfacial cell at x = %g, y = %g (t = %g).\n", x, y, t);
        dump (file = "dump-err");
        assert (0);
      }
      for (int k = 0; k < 2; k++) {
        int i = 2*k - 1;
        foreach_dimension()
          if (fs.x[k] > 0. && fs.x[k] < 1. && (cfg[i] + cfg[]) == 0.) {
            fprintf (ferr, "Inconsistent cfg distribution around x = %g, y = %g (t = %g).\n", x, y, t);
            dump (file = "dump-err");
            assert (0);
          }
      }
    }
}

void update_cfg (const scalar c, const scalar cs, face vector fs, scalar cfg)
{
  check_cfg (c, cs, fs, cfg);

  scalar cflag[];
  foreach()
    cflag[] = 0.;

  foreach()
    if ((cfg[] == 2. && c[] <= 0.) ||
        (cfg[] == -2. && c[] >= cs[])) {
        assert (cs[] > 0. && cs[] < 1.);
        cflag[] = 1.;
        for (int k = 0; k < 2; k++) {
          int i = 2*k - 1;
          foreach_dimension()
            if (fs.x[k] > 0. && fs.x[k] < 1. && cfg[i] != cfg[])
              cflag[] = 0.;
        }
    }

  foreach()
    if (cflag[] == 1.) {
      assert (cfg[] == 2. || cfg[] == -2.);
      cfg[] = - cfg[];
    }
    else if (cs[] > 0. && cs[] < 1.) { // split
      for (int k = 0; k < 2; k++) {
        int i = 2*k - 1;
        foreach_dimension()
          if (fs.x[k] > 0. && fs.x[k] < 1. && cflag[i] == 1. && cflag[] == 0.)
            cfg[] = 1.;
      }
    }

  check_cfg (c, cs, fs, cfg);

  foreach()
    if (cs[] > 0. && cs[] < 1. && cfg[] == 1.) { // merge
      int nsl = 0;
      int nsg = 0;
      for (int k = 0; k < 2; k++) {
        int i = 2*k - 1;
        foreach_dimension()
          if (fs.x[k] > 0. && fs.x[k] < 1.) {
            if (cfg[i] == 2.)
              nsl++;
            if (cfg[i] == -2.)
              nsg++;
          }
      }
      if (nsl == 2 || nsg == 2)
        cfg[] = (nsl > nsg ? 2. : -2.);
    }
}

#endif
