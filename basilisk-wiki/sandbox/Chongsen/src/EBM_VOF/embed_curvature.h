

#if TREE
static void curvature_restriction (Point point, scalar kappa)
{
  double k = 0., s = 0.;
  foreach_child()
    if (kappa[] != nodata)
      k += kappa[], s++;
  kappa[] = s ? k/s : nodata;
}



static void curvature_prolongation (Point point, scalar kappa)
{
  foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)
    #if dimension > 1
      for (int j = 0; j <= 1; j++)
    #endif
      #if dimension > 2
	for (int k = 0; k <= 1; k++)
      #endif
	  if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
	    sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    kappa[] = s ? sk/s : nodata;
  }
}
#endif // TREE


#include "embed_height_normal.h"

#include "parabola.h"

static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      foreach_dimension()
	d2 += sq(p[i].x - p[j].x);
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}

static double height_curvature_fit (Point point, scalar c, vector h)
{

  
  coord ip[dimension == 2 ? 6 : 27];
  int n = 0;


  foreach_dimension() {


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

  
  if (independents (ip, n) < (dimension == 2 ? 3 : 9))
    return nodata;

  
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
  

  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}


static double centroids_curvature_fit (Point point, scalar c)
{

  
  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (c[], m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);


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
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;
#if AXI
  parabola_fit_axi_curvature (&fit, y + fc.y*Delta, Delta, &kappa, NULL);
#endif
  return kappa;
}


typedef struct {
  int h; // number of standard HF curvatures
  int f; // number of parabolic fit HF curvatures
  int a; // number of averaged curvatures
  int c; // number of centroids fit curvatures
} cstats;

trace
cstats curvature (scalar c, scalar tmp_c, scalar mark, scalar kappa,
		  double sigma = 1.[0], bool add = false)
{
  int sh = 0, sf = 0, sa = 0, sc = 0;

  scalar interface_j[];

  foreach(){
    if (mark[]==4 || mark[]==5 || mark[]==6)
      interface_j[]=1;
    else
      interface_j[]=0;
  }

  vector ch = tmp_c.height,     h = automatic (ch);
  
  if (!ch.x.i)
    heights (c, mark, h);
  

  
#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  scalar k[];
  scalar_clone (k, kappa);

  foreach()
    
    k[]=nodata;

  
  foreach(reduction(+:sh) reduction(+:sf)) {

    if  (k[]==nodata){
    
  
    if (interface_j[]==0)
      k[] = nodata;


    else if ((k[] = height_curvature (point, tmp_c, h)) != nodata){
      sh++;
      }
      
    else if ((k[] = height_curvature_fit (point, tmp_c, h)) != nodata)
      sf++;

    }
    
  }
  
  foreach (reduction(+:sa) reduction(+:sc)) {
    

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interface_j[]==1) {

      
      double sk = 0., a = 0.;
      foreach_neighbor(1)
	if (k[] < nodata)
	  sk += k[], a++;
      if (a > 0.)
	kf = sk/a, sa++;
      else


	kf = centroids_curvature_fit (point, tmp_c), sc++;
    }
    else
      kf = nodata;

    
    if (kf == nodata)
      kappa[] = nodata;
    else if (add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;      
  }


  return (cstats){sh, sf, sa, sc};
}






/**
# Position of an interface

This is similar to curvature but this time for the position of the
interface, defined as
$$
pos = \mathbf{G}\cdot(\mathbf{x} - \mathbf{Z})
$$
with $\mathbf{G}$ and $\mathbf{Z}$ two vectors and $\mathbf{x}$ the
coordinates of the interface.

This is defined only in interfacial cells. In all the other cells it
takes the value *nodata*.

We first need a function to compute the position $\mathbf{x}$ of an
interface. For accuracy, we first try to use height functions. */

foreach_dimension()
static double pos_x (Point point, vector h, coord * G, coord * Z)
{
  if (fabs(height(h.x[])) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.x += height(h.x[])*Delta;
  double pos = 0.;
  foreach_dimension()
    pos += (o.x - Z->x)*G->x;
  return pos;
}

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the position. This is done by the function below which returns
the HF position given a volume fraction field *f*, a height function
field *h* and vectors *G* and *Z*. */

static double height_position (Point point, scalar f, vector h,
			       coord * G, coord * Z)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *f*) and corresponding HF position function *pos*
  (defined above). */

  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  foreach_dimension()
    n.x.n = f[1] - f[-1], n.x.pos = pos_x;
  
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormPos, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormPos, n.y, n.z);
#endif

  /**
  We try each position function in turn. */

  double pos = nodata;
  foreach_dimension()
    if (pos == nodata)
      pos = n.x.pos (point, h, G, Z);

  return pos;
}
