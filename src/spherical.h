/**
# Spherical coordinates

The default radius of the sphere is set to one. */

double Radius = 1.;

/**
Rather than the classical [mathematical
convention](http://en.wikipedia.org/wiki/Spherical_coordinate_system),
we use the geographic convention: *x* is the longitude within
$[-180:180]$ degrees and *y* is the latitude within $[-90:90]$
degrees. $\Delta$ is the characteristic cell length expressed in the
same units as *Radius*. */

map()
{
  x = x < -180. ? x + 360. : x > 180. ? x - 360. : x;
  Delta_x = Delta_y = Delta;
  Delta *= Radius*pi/180.;
}

/**
On trees we need to define consistent refinement functions. The
default restriction functions (averages) are already consistent
(i.e. they preserve volume and length integrals).

Mathematically we have
$$
cm = fm.y = \cos(\phi)
$$
$$
fm.x = 1
$$
with $\phi$ the latitude in radians. This is the definition we use for
the length scale factors *fm*.

For the volume scale factor *cm*, more care needs to be taken to
guarantee discrete volume conservation i.e.
$$
\sum_{children}cm_{child} = 4cm_{parent}
$$
This is achieved by computing the exact area of a cell i.e.
$$
\Delta^2 \int^{\phi + d\phi/2}_{\phi - d \phi/2} \cos\phi d\phi =
\Delta^2  (\sin(\phi + d\phi/2) - \sin(\phi - d\phi/2)) = \Delta^2 cm
$$
*/

#if TREE
static void refine_cm_lonlat (Point point, scalar cm)
{
  double phi = y*pi/180., dphi = Delta/(2.*Radius);
  double sinphi = sin(phi);
  fine(cm,0,0) = fine(cm,1,0) = (sinphi - sin(phi - dphi))/dphi;
  fine(cm,0,1) = fine(cm,1,1) = (sin(phi + dphi) - sinphi)/dphi;
}

static void refine_face_x_lonlat (Point point, scalar fm)
{
  if (!is_refined(neighbor(-1)))
    fine(fm,0,0) = fine(fm,0,1) = 1.;
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors)
    fine(fm,2,0) = fine(fm,2,1) = 1.;
  fine(fm,1,0) = fine(fm,1,1) = 1.;
}

static void refine_face_y_lonlat (Point point, scalar fm)
{
  double phi = y*pi/180., dphi = Delta/(2.*Radius);
  if (!is_refined(neighbor(0,-1)))
    fine(fm,0,0) = fine(fm,1,0) = cos(phi - dphi);
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors)
    fine(fm,0,2) = fine(fm,1,2) = cos(phi + dphi);
  fine(fm,0,1) = fine(fm,1,1) = cos(phi);
}
#endif

event metric (i = 0) {

  /**
  We initialise the scale factors, taking care to first allocate the
  fields if they are still constant. */

  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }
  scalar cmv = cm;
  foreach() {
    double phi = y*pi/180., dphi = Delta/(2.*Radius);
    cmv[] = (sin(phi + dphi) - sin(phi - dphi))/(2.*dphi);
  }

  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face(x)
    fmv.x[] = 1.;
  foreach_face(y)
    fmv.y[] = cos(y*pi/180.);

  /**
  We set our refinement/prolongation functions on trees. */

#if TREE
  cm.refine = cm.prolongation = refine_cm_lonlat;
  fm.x.prolongation = refine_face_x_lonlat;
  fm.y.prolongation = refine_face_y_lonlat;
#endif
}
