typedef double real;

#define GRIDNAME "Cartesian 1D"
#define dimension 1
#define GHOSTS 1

#define _I     (point.i - 1)
#define _DELTA (1./point.n)

typedef struct {
  Grid g;
  char * d;
  int n;
} Cartesian;

struct _Point {
  int i, j, level, n;
};
static Point last_point;

#define cartesian ((Cartesian *)grid)

@define data(k,l,m) ((double *)&cartesian->d[(point.i + k)*datasize])
@define allocated(...) true

macro POINT_VARIABLES (Point point = point) { VARIABLES(); }

macro2 foreach (char flags = 0, Reduce reductions = None)
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
    point.n = cartesian->n;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 1; _k <= point.n; _k++) {
	point.i = _k;
	{...}
      }
  }
}
  
macro2 foreach_face_generic (char flags = 0, Reduce reductions = None,
				const char * order = "xyz")
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
    point.n = cartesian->n;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 1; _k <= point.n + 1; _k++) {
	point.i = _k;
	{...}
      }
  }
}

macro1 is_face_x() {{ int ig = -1; NOT_UNUSED(ig); {...} }}

// ghost cell coordinates for each direction
static int _ig[] = {1,-1};

// Box boundaries

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  if (!list)
    return;
  
  int d = ((BoxBoundary *)b)->d;

  Point point;
  point.n = cartesian->n;
  int ig = _ig[d];
  assert (d <= left);
  point.i = d == right ? point.n + GHOSTS : GHOSTS;
  Point neighbor = {point.i + ig};
  for (scalar s in list) {
    scalar b = s.v.x;
    val(s,ig) = b.boundary[d] (point, neighbor, s, NULL);
  }
}

static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[d] != periodic_bc) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  if (centered) {
    Point point;
    point.n = cartesian->n;
    int ig = _ig[d];
    point.i = d == right ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig};
    for (scalar s in centered)
      val(s,ig) = s.boundary[d] (point, neighbor, s, NULL);
    free (centered);
  }
    
  box_boundary_level_normal (b, normal, l);
  free (normal);
}

// periodic boundaries

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  Point point = *((Point *)grid);
  point.i = 0, point.n = N;
  for (int i = 0; i < GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i + point.n];
  for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i - point.n];

  free (list1);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  free (cartesian->d);
  free (cartesian);
  grid = NULL;
}

@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  char * data = cartesian->d;
  for (int i = 0; i < cartesian->n + 2; i++, data += datasize) {
    double * v = (double *) data;
    for (scalar s in list)
      if (!is_constant(s))
	v[s.i] = val;
  }
}

void init_grid (int n)
{
  if (cartesian && n == cartesian->n)
    return;
  free_grid();
  Cartesian * p = qmalloc (1, Cartesian);
  size_t len = (n + 2)*datasize;
  p->n = N = n;
  p->d = qmalloc (len, char);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->d;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = (Grid *) p;
  reset (all, 0.);
  // box boundaries
  for (int d = 0; d < 2; d++) {
    BoxBoundary * box = qcalloc (1, BoxBoundary);
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
    add_boundary (b);
  }
  // periodic boundaries
  Boundary * b = qcalloc (1, Boundary);
  b->level = periodic_boundary_level_x;
  add_boundary (b);
  // mesh size
  grid->n = grid->tn = n;
}

void realloc_scalar (int size)
{
  Cartesian * p = cartesian;
  size_t len = (p->n + 2);
  qrealloc (p->d, len*(datasize + size), char);
  char * data = p->d + (len - 1)*datasize;
  for (int i = p->n + 1; i > 0; i--, data -= datasize)
    memmove (data + i*size, data, datasize);
  datasize += size;
}

Point locate (double xp = 0, double yp = 0, double zp = 0)
{
  Point point;
  point.n = cartesian->n;
  double a = (xp - X0)/L0*point.n;
  point.i = a + 1;
  point.level = (a > -0.5 && a < point.n + 0.5) ? 0 : - 1;
  return point;
}

#include "variables.h"
#include "cartesian-common.h"

macro2 foreach_vertex (char flags = 0, Reduce reductions = None)
{
  foreach_face_generic (flags, reductions) {
    int ig = -1; NOT_UNUSED(ig);
    {...}
  }
}

void cartesian1D_methods()
{
  cartesian_methods();
}
