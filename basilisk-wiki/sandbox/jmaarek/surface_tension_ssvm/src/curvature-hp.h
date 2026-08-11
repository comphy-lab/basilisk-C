/**
# Curvature of an interface

The curvature field is defined only in interfacial cells. In all the
other cells it takes the value *nodata*. 

On trees, we need to redefine the restriction function to take
this into account i.e. the curvature of the parent cell is the average
of the curvatures in the interfacial child cells. */


#if dimension > 1

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
orientation. */


#if dimension == 2
foreach_dimension()
static double kappa_hp_y (Point point, vector h)
{
  int ori = orientation(h.y[]);
  for (int i = -2; i <= 2; i++)
    if (h.y[i] == nodata || orientation(h.y[i]) != ori)
      return nodata;
  double weight = 0.0; //complete suppression of Nyquist mode in longitunal direction
  double hx = (h.y[1] - h.y[-1])/2.;
  double hxx = (1-weight)*(h.y[-2] - 2.*h.y[] + h.y[2])/(4.*Delta)
           + weight*(h.y[-1] - 2.*h.y[] + h.y[1])/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}

#else // dimension == 3
foreach_dimension()
static double kappa_hp_z (Point point, vector h)
{
  int ori = orientation(h.z[]);

  // Require full 5x5 stencil — nodata on any point returns nodata cleanly
  for (int i = -2; i <= 2; i++)
    for (int j = -2; j <= 2; j++)
      if (h.z[i,j] == nodata || orientation(h.z[i,j]) != ori)
        return nodata;

  // First derivatives: use 2*Delta spacing for consistency with hxx, hyy
  double hx = (h.z[1] - h.z[-1])/2.;
  double hy = (h.z[0,1] - h.z[0,-1])/2.;

  // Second derivatives: longitudinal stencil uses h[±2,j] to eliminate
  // Nyquist in the derivative direction. Transverse averaging over
  // j = -1, 0, 1 retained from original formulation.
  double filter = 0.5; //switched from 0.2 //complete suppression of Nyquist mode in transverse direction
  double hxx = (filter*(h.z[ 2, 1] + h.z[-2, 1] - 2.*h.z[0, 1]) +
                       (h.z[ 2   ] + h.z[-2   ] - 2.*h.z[     ]) +
                filter*(h.z[ 2,-1] + h.z[-2,-1] - 2.*h.z[0,-1]))/
    ((1. + 2.*filter)*4.*Delta);
  double hyy = (filter*(h.z[ 1, 2] + h.z[ 1,-2] - 2.*h.z[ 1]) +
                       (h.z[ 0, 2] + h.z[ 0,-2] - 2.*h.z[  ]) +
                filter*(h.z[-1, 2] + h.z[-1,-2] - 2.*h.z[-1]))/
    ((1. + 2.*filter)*4.*Delta);

  // Cross derivative: consistent 2*Delta spacing in both directions
  //double hxy = (h.z[2,2] + h.z[-2,-2] - h.z[2,-2] - h.z[-2,2])/(16.*Delta);
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}

#endif

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */

static double height_curvature_hp (Point point, scalar c, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding HF curvature function *kappa*
  (defined above). */

  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;  
  // we use the smooth version to compute curvatures.
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = kappa_hp_x;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#if dimension == 3
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
	if (n.x.n < 0.)
	  kappa = - kappa;
      }
    }

  if (kappa != nodata) {
    
    /**
    We limit the maximum curvature to $0.5/\Delta$.
    The threshold has been changed compared with the heigh_curvature function. */
    if (fabs(kappa) > 0.5/Delta)
      kappa = 0.5*sign(kappa)/Delta;
    
    /**
    We add the axisymmetric curvature if necessary. */
      
#if AXI
    double nr, r = y, hx;
    if (kappaf == kappa_hp_x) {
      hx = (height(h.x[0,1]) - height(h.x[0,-1]))/2.;
      nr = hx*(orientation(h.x[]) ? 1 : -1);
    }
    else {
      r += height(h.y[])*Delta;
      hx = (height(h.y[1,0]) - height(h.y[-1,0]))/2.;
      nr = orientation(h.y[]) ? -1 : 1;
    }
    /* limit the minimum radius to half the grid size */
    kappa += nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
#endif
  }
  
  return kappa;
}


/**
Given a volume fraction field *c* and a height function field *h*,
this function returns the "mixed heights" parabola-fitted curvature
(or *nodata* if the curvature cannot be computed). */

static double height_curvature_fit_hp (Point point, scalar c, vector h)
{

  /**
  The coordinates of the interface points and the number of
  interface points. */

  coord ip[dimension == 2 ? 10 : 75];
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
    for (int i = -2; i <= 2; i++)
      if (h.y[i] != nodata && orientation(h.y[i]) == ori)
	ip[n].x = i, ip[n++].y = height(h.y[i]);
#else // dimension == 3
    for (int i = -2; i <= 2; i++)
      for (int j = -2; j <= 2; j++)
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
  double kappa = parabola_fit_curvature (&fit, 0.5, NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

/**
## Parabolic fit of centroids

If all else fails, we try a parabolic fit of interface centroids. */

static double centroids_curvature_fit_hp (Point point, scalar c)
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
  double kappa = parabola_fit_curvature (&fit, 0.5, NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

#endif // dimension > 1

/**
## General curvature computation

We first need to define "interfacial cells" i.e. cells which contain
an interface. A simple test would just be that the volume fraction is
neither zero nor one. As usual things are more complicated because of
round-off errors. They can cause the interface to be exactly aligned
with cell boundaries, so that cells on either side of this interface
have fractions exactly equal to zero or one. The function below takes
this into account. */


/**
The function below computes the mean curvature *kappa* of the
interface defined by the volume fraction *c*. It uses a combination of
the methods above: statistics on the number of curvatures computed
which each method is returned in a *cstats* data structure. 

The curvature is multiplied by *sigma* (default is one).

If *add* is *true*, the curvature is added to field *kappa*. */

trace
cstats curvature_hp (scalar c, scalar kappa,
		  double sigma = 1.[0], bool add = false)
{
  int sh = 0, sf = 0, sa = 0, sc = 0;

  /**
  On trees we set the prolongation and restriction functions for
  the curvature. */
  
#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

#if dimension > 1
  
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);

  foreach(reduction(+:sh) reduction(+:sf)) {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */

    if (!interfacial (point, c))
      k[] = nodata;

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_curvature_hp (point, c, h)) != nodata)
      sh++;
    else if ((k[] = height_curvature_fit_hp (point, c, h)) != nodata)
      sf++;
  }
  
  foreach (reduction(+:sa) reduction(+:sc)) {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial (point, c)) {

      /**
      ...or the average of the curvatures in the $3^{d}$ neighborhood
      of interfacial cells. */
      
      double sk = 0., a = 0.;
      foreach_neighbor(1)
	if (k[] < nodata)
	  sk += k[], a++;
      if (a > 0.)
	kf = sk/a, sa++;
      else

	/**
	Empty neighborhood: we try centroids as a last resort. */

	kf = centroids_curvature_fit_hp (point, c), sc++;
    }
    else
      kf = nodata;

    /**
    We add or set *kappa*. */
    
    if (kf == nodata)
      kappa[] = nodata;
    else if (add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;      
  }
  
#else // dimension == 1
  foreach() {
    if (!interfacial (point, c))
      kappa[] = nodata;
    else {
      double r = x + sign(c[-1] - c[1])*(clamp(c[],0.,1.) - 0.5)*Delta;
      double p = r > 0. ? - 2.*sigma/r : 0.;
      if (add)
	kappa[] += p;
      else
	kappa[] = p;
    }
  }
#endif // dimension == 1

  return (cstats){sh, sf, sa, sc};   
}
