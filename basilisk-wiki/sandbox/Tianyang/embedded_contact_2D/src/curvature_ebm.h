#include "heights_contact_ebm.h"

#if dimension > 1

coord interface_normal (Point point, scalar c);

static double height_curvature_fit_ebm (Point point, scalar c, vector h)
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

  scalar cet = c.cet;
  assert (cet[] != nodata);

  coord m = interface_normal (point, c), fc;
  double alpha = plane_alpha (cet[], m);
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
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

static double centroids_curvature_fit_ebm (Point point, scalar c)
{

  /**
  We recover the interface normal and the centroid of the interface
  fragment and initialize the parabolic fit. */

  scalar cet = c.cet;
  assert (cet[] != nodata);

  coord m = interface_normal (point, c), fc;
  double alpha = plane_alpha (cet[], m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  /**
  We add the interface centroids in a $3^d$ neighborhood and compute
  the curvature. */

  coord r = {x,y,z};
  foreach_neighbor(1)
    if (cet[] > 0. && cet[] < 1.) {
      coord m = interface_normal (point, c), fc;
      double alpha = plane_alpha (cet[], m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      foreach_dimension()
	fc.x += (rn.x - r.x)/Delta;
      parabola_fit_add (&fit, fc, area);
    }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}

#endif // dimension > 1

static inline bool interfacial_ebm (Point point, scalar c, scalar cs)
{
  if (cs[] <= 0.)
    return false;

  if (c[] >= cs[]) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] <= 0. && cs[i] > 0.)
	  return true;
  }
  else if (c[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] >= cs[i] && cs[i] > 0.)
	  return true;
  }
  else // c[] > 0. && c[] < cs[]
    return true;
  return false;
}

trace
cstats curvature_ebm (scalar c, scalar cs, scalar kappa,
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
  assert (ch.x.i);

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);

  foreach(reduction(+:sh) reduction(+:sf)) {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */

    if (!interfacial_ebm (point, c, cs))
      k[] = nodata;

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_curvature (point, c, h)) != nodata)
      sh++;
    else if ((k[] = height_curvature_fit_ebm (point, c, h)) != nodata)
      sf++;
  }
  
  foreach (reduction(+:sa) reduction(+:sc)) {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial_ebm (point, c, cs)) {

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

	kf = centroids_curvature_fit_ebm (point, c), sc++;
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

#endif // dimension == 1

  if (sf + sa + sc > 0)
    fprintf (stdout, "Warning: height_curvature fails (sf = %d, sa = %d, sc = %d) at t = %g!\n", sf, sa, sc, t);

  return (cstats){sh, sf, sa, sc};   
}
