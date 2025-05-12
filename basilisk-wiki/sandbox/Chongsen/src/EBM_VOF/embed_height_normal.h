#include "embed_heights.h"


#if dimension == 3
foreach_dimension()
static double kappa_z (Point point, vector h)
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
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}


foreach_dimension()
static double kappa1_z (Point point, vector h)
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
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}

foreach_dimension()
static coord normal2_z (Point point, vector h)
{
  scalar hz = h.z;
  if (hz[] == nodata)
    return (coord){nodata, nodata, nodata};
  int ori = orientation(hz[]);
  double a = ori ? -1. : 1.;
  coord n;
  n.z = a;
  foreach_dimension(2) {
    if (allocated(-1) && hz[-1] != nodata && orientation(hz[-1]) == ori) {
      if (allocated(1) && hz[1] != nodata && orientation(hz[1]) == ori)
	n.x = a*(hz[-1] - hz[1])/2.;
      else
	n.x = a*(hz[-1] - hz[]);
    }
    else if (allocated(1) && hz[1] != nodata && orientation(hz[1]) == ori)
      n.x = a*(hz[] - hz[1]);
    else
      n.x = nodata;
  }
  return n;
}

foreach_dimension()
static coord normal_z (Point point, vector h) {
  coord n = normal2_z (point, h);
  double nn = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (nn < nodata) {
    foreach_dimension()
      n.x /= nn;
    return n;
  }
  return (coord){nodata, nodata, nodata};
}


#else // dimension == 2
foreach_dimension()
static double kappa_y (Point point, vector h)
{
  int ori = orientation(h.y[]);
  for (int i = -1; i <= 1; i++)

    if (h.y[i] == nodata || orientation(h.y[i]) != ori)
      return nodata;
  double hx = (h.y[1] - h.y[-1])/2.;
  double hxx = (h.y[1] + h.y[-1] - 2.*h.y[])/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}

foreach_dimension()
static coord normal_y (Point point, vector h)
{
  coord n = {nodata, nodata, nodata};
  if (h.y[] == nodata)
    return n;
  int ori = orientation(h.y[]);


    if (h.y[-1] != nodata && orientation(h.y[-1]) == ori) {
    if (h.y[1] != nodata && orientation(h.y[1]) == ori)
      n.x = (h.y[-1] - h.y[1])/2.;
    else
      n.x = h.y[-1] - h.y[];
  }
  else if (h.y[1] != nodata && orientation(h.y[1]) == ori)
    n.x = h.y[] - h.y[1];
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;
  
}



#endif

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */

static double height_curvature (Point point, scalar c, vector h)
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
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = kappa_x;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  coord ns = {0.5,0.5};
  ns.x = fs.x[]-fs.x[1];
  ns.y = fs.y[]-fs.y[0,1];
  if (ns.x==0&&ns.y==0){
    ns.x = 0.5;
    ns.y = 0.5;
  }

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
     We limit the maximum curvature to $1/\Delta$. */
	
    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
    
    /**
     We add the axisymmetric curvature if necessary. */

#if AXI
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
    kappa += nr/max (sqrt(1. + sq(hx))*r, Delta/2.);

#endif
  }
  
  return kappa;
}


/**
The function below works in a similar manner to return the normal
estimated using height-functions (or a *nodata* vector if this cannot
be done). */
coord height_normal (Point point, scalar c, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding normal function *normal*
  (defined above). */

  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;  
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.normal = normal_x;
  
  /**
  We sort these pairs in decreasing order of $|n|$. */

  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormNormal, n.x, n.y);

#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormNormal, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormNormal, n.y, n.z);
#endif

  /**
  We try each normal function in turn. */

  coord normal = {nodata, nodata, nodata};
  foreach_dimension()
    if (normal.x == nodata)
      normal = n.x.normal (point, h);
  
  return normal;
}


/**
In three dimensions, these functions return the (two) components of
the normal projected onto the $(x,y)$ plane (respectively). */

#if dimension == 3
foreach_dimension()
coord height_normal_z (Point point, vector h)
{
  coord nx = normal2_x (point, h);
  coord ny = normal2_y (point, h);
  if (fabs(nx.y) < fabs(ny.x)) {
    normalize (&nx);
    return nx;
  }
  else if (ny.x != nodata) {
    normalize (&ny);
    return ny;
  }
  return (coord){nodata, nodata, nodata};
}
#endif



extern scalar * interfaces_mark;


coord interface_normal (Point point, scalar tmp_c);
#define interface_normal(point, tmp_c) interface_normal (point, tmp_c)


#undef SEPS
#define SEPS 1e-30

#ifndef VFTL
# define VFTL 1e-10
#endif

#include "fractions.h"

coord interface_normal (Point point, scalar tmp_c)
{
  coord n = {nodata, nodata, nodata};
  



  if (!tmp_c.height.x.i || (n = height_normal (point, tmp_c, tmp_c.height)).x == nodata){
    n = mycs (point, tmp_c);
  }
  

  return n;
}





