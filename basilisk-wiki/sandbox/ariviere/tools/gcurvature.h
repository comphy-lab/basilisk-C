/**
Small modification of curvature.h of the main branch in order to compute the gaussian curvature. 

# Curvature of an interface
Add the gaussian curvature in 3d only. In 2D it should be 0 but not tested in this case.

The curvature field is defined only in interfacial cells. In all the
other cells it takes the value *nodata*. 

On trees, we need to redefine the restriction function to take
this into account i.e. the curvature of the parent cell is the average
of the curvatures in the interfacial child cells. */


/**
## Height-function curvature and normal

To compute the curvature, we estimate the derivatives of the height
functions in a given direction (*x*, *y* or *z*). We first check that
all the heights are defined and that their orientations are the
same. We then compute the curvature as
$$
\kappa = \frac{h_{xx}}{(1 + h_x^2)^{3/2}}
$$
in two dimensions, or
$$
\kappa = \frac{h_{xx}(1 + h_y^2) + h_{yy}(1 + h_x^2) - 2h_{xy}h_xh_y}
{(1 + h_x^2 + h_y^2)^{3/2}}
$$
in three dimensions.

The normal is computed in a similar way, but also allowing for
asymmetric 2-points stencils and taking into account the
orientation. 
 
The Gaussian curvature in 3D is computed with the following formula:
$$ 
G = \frac{h_{xx} h_{yy} - h_{xy}^2}{(1+h_x^2 + h_y^2)^2}
$$
 
 */

#include "heights.h"

#if dimension==3
foreach_dimension()
static double gkappa_z (Point point, vector h)
{
  int ori = orientation(h.z[]);
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (h.z[i,j] == nodata || orientation(h.z[i,j]) != ori)
	return nodata;
  double hx = (h.z[1] - h.z[-1])/2.;
  double hy = (h.z[0,1] - h.z[0,-1])/2.;

  /**
  We "filter" the curvature using a weighted sum of the three
  second-derivatives in the $x$ and $y$ directions. This is necessary
  to avoid a numerical mode when the curvature is used to compute
  surface tension. */
  
  double filter = 0.2;
  double hxx = (filter*(h.z[1,1] + h.z[-1,1] - 2.*h.z[0,1]) +
		(h.z[1] + h.z[-1] - 2.*h.z[]) +
		filter*(h.z[1,-1] + h.z[-1,-1] - 2.*h.z[0,-1]))/
    ((1. + 2.*filter)*Delta);
  double hyy = (filter*(h.z[1,1] + h.z[1,-1] - 2.*h.z[1]) +
		(h.z[0,1] + h.z[0,-1] - 2.*h.z[]) +
		filter*(h.z[-1,1] + h.z[-1,-1] - 2.*h.z[-1]))/
    ((1. + 2.*filter)*Delta);
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*Delta);
  return (hxx*hyy - sq(hxy))/sq(1+sq(hx)+sq(hy));
}

#endif

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */


static double height_gcurvature (Point point, scalar c, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding HF gcurvature function *kappa*
  (defined above). */

  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;  
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = gkappa_x;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#if dimension==3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);
#endif

  /**
  We try each curvature function in turn. */

  double kappa = nodata;
  foreach_dimension()
    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
	kappaf = n.x.kappa;
//	if (n.x.n < 0.) //sign must NOT change when change the norm!
//	  kappa = - kappa;
      }
    }

  if (kappa != nodata) {
    
    /**
     We limit the maximum gaussian curvature to $(1/\Delta)^2. (Remember the
     relation between mean and gaussian curvature : $\kappa^2 \geq \kappa_g$ $. */
	
    if (fabs(kappa) > 1./sq(Delta))
      kappa = sign(kappa)/sq(Delta);
    
    /**
     We add the axisymmetric curvature if necessary. */
      
#if AXI //not sure but I just multiplied by the new curvature
    double nr, r = y, hx;
    if (kappaf == kappa_x) {
      hx = (height(h.x[0,1]) - height(h.x[0,-1]))/2.;
      nr = hx*(orientation(h.x[]) ? 1 : -1);
    }
    else {
      r += height(h.y[])*Delta;
      hx = (height(h.y[1,0]) - height(h.y[-1,0]))/2.;
      nr = orientation(h.y[]) ? -1 : 1;
    }
    /* limit the minimum radius to half the grid size */
    kappa *= nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
#endif
  }
  
  return kappa;
}

/**
The function below works in a similar manner to return the normal
estimated using height-functions (or a *nodata* vector if this cannot
be done). */


/**
## Parabolic fit of "mixed" height-functions

When the standard height function curvature calculation is not
possible (for example because not enough heights are available in any
given direction), one can try to combine all the available heights
(thus using "mixed" directions) to obtain points on the
interface. These point locations can then be fitted with a parabola
(using least-mean-square optimisation) and the resulting curvature can
be computed. The fitting functions are defined in the file included
below. */


/**
Given a volume fraction field *c* and a height function field *h*,
this function returns the "mixed heights" parabola-fitted curvature
(or *nodata* if the curvature cannot be computed). */


#include "gparabola.h"

static double height_gcurvature_fit (Point point, scalar c, vector h)
{

  /**
  The coordinates of the interface points and the number of
  interface points. */
  
  coord ip[dimension == 2 ? 6 : 27];
  int n = 0;

  /**
  We collect the points along all directions. */
  
  foreach_dimension() {

    /**
    We don't want to mix heights with different orientations. We first
    find the "dominant" orientation *ori*. */
    
    int n1 = 0, n2 = 0;
#if dimension == 2
    for (int i = -1; i <= 1; i++)
      if (h.y[i] != nodata) {
	if (orientation(h.y[i])) n1++; else n2++;
      }
#else // dimension == 3
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h.z[i,j] != nodata) {
	  if (orientation(h.z[i,j])) n1++; else n2++;
	}
#endif
    int ori = (n1 > n2);

    /**
    We look for height-functions with the dominant orientation and
    store the corresponding interface coordinates (relative to the
    center of the cell and normalised by the cell size). */

#if dimension == 2
    for (int i = -1; i <= 1; i++)
      if (h.y[i] != nodata && orientation(h.y[i]) == ori)
	ip[n].x = i, ip[n++].y = height(h.y[i]);
#else // dimension == 3
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h.z[i,j] != nodata && orientation(h.z[i,j]) == ori)
	  ip[n].x = i, ip[n].y = j, ip[n++].z = height(h.z[i,j]);
#endif
  }

  /**
  If we don't have enough independent points, we cannot do the
  parabolic fit. */
  
  if (independents (ip, n) < (dimension == 2 ? 3 : 9))
    return nodata;

  /**
  We recover the interface normal and the centroid of the interface
  fragment and initialize the parabolic fit. */
  
  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (c[], m);
  double area = plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);
#if dimension == 2
  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, PARABOLA_FIT_CENTER_WEIGHT);
#else // dimension == 3
  parabola_fit_add (&fit, fc, area*100.);
#endif
  
  /**
  We add the collected interface positions and compute the
  curvature. */

  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_gcurvature (&fit, sq(2.))/sq(Delta);
#if AXI//TO CHECK
  parabola_fit_axi_gcurvature (&fit, y + fc.y*Delta, Delta, &kappa);
#endif
  return kappa;
}

/**
## Parabolic fit of centroids

If all else fails, we try a parabolic fit of interface centroids. */

static double centroids_gcurvature_fit (Point point, scalar c)
{

  /**
  We recover the interface normal and the centroid of the interface
  fragment and initialize the parabolic fit. */
  
  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (c[], m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  /**
  We add the interface centroids in a $3^d$ neighborhood and compute
  the curvature. */

  coord r = {x,y,z};
  foreach_neighbor(1)
    if (c[] > 0. && c[] < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = plane_alpha (c[], m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      foreach_dimension()
	fc.x += (rn.x - r.x)/Delta;
      parabola_fit_add (&fit, fc, area);
    }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_gcurvature (&fit, sq(2.))/sq(Delta);
#if AXI
  parabola_fit_axi_gcurvature (&fit, y + fc.y*Delta, Delta, &kappa);
#endif
  return kappa;
}

/**
## General curvature computation

The function below computes the gaussian curvature *kappa* of the
interface defined by the volume fraction *c*. It uses a combination of
the methods above: statistics on the number of curvatures computed
which each method is returned in a *cstats* data structure. 

