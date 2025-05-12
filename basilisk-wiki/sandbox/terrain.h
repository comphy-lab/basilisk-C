/**
This uses the kdt module to reconstruct the terrain from the list of kdt databases specified
*/
#include <stdarg.h>
#include <kdt/kdt.h>
#ifdef _OPENMP
# define NPROC omp_get_max_threads()
#else
# define NPROC 1
#endif

typedef struct {
  Kdt *** kdt;
  scalar n, dmin, dmax;
} Terrain;

Terrain * _terrain = NULL;
int _nterrain = 0;

/** In order to use the kdt_query_sum it is necessary to specify two functions, an includes function and an intersects function.
  
The includes function determines whether the cell specified by point p (with centre (x,y) and cell size delta includes a given kdt rectangle
*/ 
static int includes (KdtRect rect, Point * p)
{
  Point point = *p;
  delta /= 2.;
  return (rect[0].l >= x - delta && rect[0].h <= x + delta &&
	  rect[1].l >= y - delta && rect[1].h <= y + delta);
}
/**
The intersects function determines whether there is any intersection between the cell specified by point p and the kdt rectangle
*/
static int intersects (KdtRect rect, Point * p)
{
  Point point = *p;
  delta /= 2.;
  return (rect[0].l <= x + delta && rect[0].h >= x - delta &&
	  rect[1].l <= y + delta && rect[1].h >= y - delta);
}

static void reconstruct_terrain (Point point, scalar zb)
{
  KdtSum s;
  kdt_sum_init (&s);
  delta /= 2.;
  KdtRect rect = {{x - delta, x + delta},
		  {y - delta, y + delta}};
  for (Kdt ** kdt = _terrain[zb].kdt[pid()]; *kdt; kdt++)
    kdt_query_sum (*kdt,
		   (KdtCheck) includes,
		   (KdtCheck) intersects, &point,
		   rect, &s);
  val(_terrain[zb].n,0,0) = s.n;
  if (s.w > 0.) {
    zb[] = s.H0/s.w;
    val(_terrain[zb].dmin,0,0) = s.Hmin;
    val(_terrain[zb].dmax,0,0) = s.Hmax;
  }
  else {
    /* not enough points in database, use bilinear interpolation
       from coarser level instead */
    if (level > 0)
      zb[] = (9.*coarse(zb,0,0) + 
	      3.*(coarse(zb,child.x,0) + coarse(zb,0,child.y)) + 
	      coarse(zb,child.x,child.y))/16.;
    else
      zb[] = 0.; // no points at level 0!
    val(_terrain[zb].dmin,0,0) = nodata;
    val(_terrain[zb].dmax,0,0) = nodata;
  }
}

static void refine_terrain (Point point, scalar zb)
{
  foreach_child()
    reconstruct_terrain (point, zb);
}

void terrain (scalar zb, ...)
{  
  if (zb >= _nterrain) {
    _terrain = realloc (_terrain, sizeof (Terrain)*(zb + 1));
    _nterrain = zb + 1;
  }
  _terrain[zb].kdt = calloc (NPROC, sizeof (Kdt **));

  int nt = 0;
  va_list ap;
  va_start (ap, zb);
  char * name;
  while ((name = va_arg (ap, char *)))
    for (int i = 0; i < NPROC; i++) {
      Kdt ** kdt = _terrain[zb].kdt[i];
      _terrain[zb].kdt[i] = kdt = realloc (kdt, sizeof(Kdt *)*(nt + 2));
      kdt[nt] = kdt_new();
      kdt[nt + 1] = NULL;
      if (kdt_open (kdt[nt], name)) {
	fprintf (stderr, "terrain: could not open terrain database '%s'\n", 
		 name);
	exit (1);
      }
    }
  va_end (ap);

  zb.refine = refine_terrain;
  scalar n = new scalar;
  n.refine = refine_none;
  n.coarsen = refine_none;
  _terrain[zb].n = n;
  scalar dmin = new scalar;
  dmin.refine = refine_none;
  dmin.coarsen = refine_none;
  _terrain[zb].dmin = dmin;
  scalar dmax = new scalar;
  dmax.refine = refine_none;
  dmax.coarsen = refine_none;
  _terrain[zb].dmax = dmax;

  trash ({zb, n, dmin, dmax});
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      reconstruct_terrain (point, zb);
    boundary_level ({zb}, l);
  }
}
