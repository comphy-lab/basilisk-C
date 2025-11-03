#if SINGLE_PRECISION
typedef float real;
#else
typedef double real;
#endif

#ifndef GRIDNAME
# define GRIDNAME "Multigrid"
#endif
#define GHOSTS 2

/* By default only one layer of ghost cells is used on the boundary to
   optimise the cost of boundary conditions. */

#ifndef BGHOSTS
@ define BGHOSTS 1
#endif

#if _MPI
# define ND(i) ((size_t)(1 << point.level))
#else
# define ND(i) ((size_t)(1 << point.level)*((int *)&Dimensions)[i])
#endif

#define _I     (point.i - GHOSTS)
#define _J     (point.j - GHOSTS)
#define _K     (point.k - GHOSTS)

int Dimensions_scale = 1;
#define _DELTA (1./((1 << point.level)*Dimensions_scale))

typedef struct {
  Grid g;
  char * d;
  size_t * shift;
} Multigrid;

struct _Point {
  int i;
#if dimension > 1
  int j;
#endif
#if dimension > 2
  int k;
#endif
  int level;
#if dimension == 1
  struct { int x; } n;
#elif dimension == 2
  struct { int x, y; } n;
#elif dimension == 3
  struct { int x, y, z; } n;
#endif
#if LAYERS
  int l;
  @define _BLOCK_INDEX , point.l
#else
  @define _BLOCK_INDEX
#endif
};
static Point last_point;

#if LAYERS
# include "grid/layers.h"
#endif

#define multigrid ((Multigrid *)grid)
#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
#if dimension == 1
@undef val
@def val(a,k,l,m) (((real *)multigrid->d)[point.i + (k) +
					  multigrid->shift[point.level] +
					  _index(a,m)*multigrid->shift[depth() + 1]])
@
#elif dimension == 2
@undef val
@def val(a,k,l,m) (((real *)multigrid->d)[point.j + (l) +
					  (point.i + (k))*(ND(1) + 2*GHOSTS) +
					  multigrid->shift[point.level] +
					  _index(a,m)*multigrid->shift[depth() + 1]])
@
#elif dimension == 3
@undef val
@def val(a,l,m,o) (((real *)multigrid->d)[point.k + (o) +
					  (ND(2) + 2*GHOSTS)*
					  (point.j + (m) +
					   (point.i + (l))*(ND(1) + 2*GHOSTS)) +
					  multigrid->shift[point.level] +
					  _index(a,0)*multigrid->shift[depth() + 1]])
@
#endif

/* low-level memory management */
#if dimension == 1
# if BGHOSTS == 1
@define allocated(...) true
# else // BGHOST != 1
@define allocated(k,l,m) (point.i+(k) >= 0 && point.i+(k) < ND(0) + 2*GHOSTS)
# endif // BGHOST != 1
@def allocated_child(k,l,m) (level < depth() &&
                             point.i > 0 && point.i <= ND(0) + 2)
@
#elif dimension == 2
# if BGHOSTS == 1
@define allocated(...) true
# else // BGHOST != 1
@def allocated(k,l,m) (point.i+(k) >= 0 && point.i+(k) < ND(0) + 2*GHOSTS &&
		       point.j+(l) >= 0 && point.j+(l) < ND(1) + 2*GHOSTS)
@
# endif // BGHOST != 1
@def allocated_child(k,l,m)  (level < depth() &&
			      point.i > 0 && point.i <= ND(0) + 2 &&
			      point.j > 0 && point.j <= ND(1) + 2)
@			   
#else // dimension == 3
# if BGHOSTS == 1
@define allocated(...) true
#else // BGHOST != 1
@def allocated(a,l,m) (point.i+(a) >= 0 &&
		       point.i+(a) < ND(0) + 2*GHOSTS &&
		       point.j+(l) >= 0 &&
		       point.j+(l) < ND(1) + 2*GHOSTS &&
		       point.k+(m) >= 0 &&
		       point.k+(m) < ND(2) + 2*GHOSTS)
@
#endif // BGHOST != 1
@def allocated_child(a,l,m)  (level < depth() &&
			      point.i > 0 && point.i <= ND(0) + 2 &&
			      point.j > 0 && point.j <= ND(1) + 2 &&
			      point.k > 0 && point.k <= ND(2) + 2)
@
#endif // dimension == 3

/***** Multigrid variables and macros *****/
@define depth()       (grid->depth)
#if dimension == 1
@def fine(a,k,l,m)
(((real *)multigrid->d)[2*point.i - GHOSTS + (k) +
			multigrid->shift[point.level + 1] +
			_index(a,m)*multigrid->shift[depth() + 1]])
@
@def coarse(a,k,l,m)
(((real *)multigrid->d)[(point.i + GHOSTS)/2 + (k) +
			multigrid->shift[point.level - 1] +
			_index(a,m)*multigrid->shift[depth() + 1]])
@

macro POINT_VARIABLES (Point point = point)
{
  VARIABLES();
  int level = point.level; NOT_UNUSED(level);
  struct { int x; } child = { 2*((point.i+GHOSTS)%2)-1 }; NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
}

#elif dimension == 2
@def fine(a,k,l,m)
(((real *)multigrid->d)[2*point.j - GHOSTS + (l) +
			(2*point.i - GHOSTS + (k))*(ND(1)*2 + 2*GHOSTS) +
			multigrid->shift[point.level + 1] +
			_index(a,m)*multigrid->shift[depth() + 1]])
@
@def coarse(a,k,l,m)
(((real *)multigrid->d)[(point.j + GHOSTS)/2 + (l) +
			((point.i + GHOSTS)/2 + (k))*(ND(1)/2 + 2*GHOSTS) +
			multigrid->shift[point.level - 1] +
			_index(a,m)*multigrid->shift[depth() + 1]])
@

macro POINT_VARIABLES (Point point = point) {
  VARIABLES();
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2; parent.j = (point.j + GHOSTS)/2;
}

#elif dimension == 3
@def fine(a,l,m,o)
(((real *)multigrid->d)[2*point.k - GHOSTS + (o) +
			(ND(2)*2 + 2*GHOSTS)*
			(2*point.j - GHOSTS + (m) +
			 (2*point.i - GHOSTS + (l))*(ND(1)*2 + 2*GHOSTS)) +
			multigrid->shift[point.level + 1] +
			_index(a,0)*multigrid->shift[depth() + 1]])
@
@def coarse(a,l,m,o)
(((real *)multigrid->d)[(point.k + GHOSTS)/2 + (o) +
			(ND(2)/2 + 2*GHOSTS)*
			((point.j + GHOSTS)/2 + (m) +
			 ((point.i + GHOSTS)/2 + (l))*(ND(1)/2 + 2*GHOSTS)) +
			multigrid->shift[point.level - 1] +
			_index(a,0)*multigrid->shift[depth() + 1]])
@

macro POINT_VARIABLES (Point point = point)
{
  VARIABLES();
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y, z; } child = {
    2*((point.i + GHOSTS)%2) - 1,
    2*((point.j + GHOSTS)%2) - 1,
    2*((point.k + GHOSTS)%2) - 1
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
  parent.j = (point.j + GHOSTS)/2;
  parent.k = (point.k + GHOSTS)/2;
}

#endif // dimension == 3

#if _MPI
#if dimension == 1
# define SET_DIMENSIONS() point.n.x = 1 << point.level
#elif dimension == 2
# define SET_DIMENSIONS() point.n.x = point.n.y = 1 << point.level
#elif dimension == 3
# define SET_DIMENSIONS() point.n.x = point.n.y = point.n.z = 1 << point.level
#endif
#else // !_MPI
#if dimension == 1
# define SET_DIMENSIONS() point.n.x = (1 << point.level)*Dimensions.x
#elif dimension == 2
# define SET_DIMENSIONS()		       \
  point.n.x = (1 << point.level)*Dimensions.x, \
  point.n.y = (1 << point.level)*Dimensions.y
#elif dimension == 3
# define SET_DIMENSIONS()		       \
  point.n.x = (1 << point.level)*Dimensions.x, \
  point.n.y = (1 << point.level)*Dimensions.y, \
  point.n.z = (1 << point.level)*Dimensions.z
#endif  
#endif // !_MPI

macro2 foreach_level (int l, char flags = 0, Reduce reductions = None) {
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    SET_DIMENSIONS();
    int _k;
    OMP(omp for schedule(static))
    for (_k = GHOSTS; _k < point.n.x + GHOSTS; _k++) {
      point.i = _k;
#if dimension > 1
      for (point.j = GHOSTS; point.j < point.n.y + GHOSTS; point.j++)
#if dimension > 2
	for (point.k = GHOSTS; point.k < point.n.z + GHOSTS; point.k++)
#endif
#endif
	  {...}
    }
  }
}

macro2 foreach (char flags = 0, Reduce reductions = None) {
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    OMP(omp for schedule(static))
      for (_k = GHOSTS; _k < point.n.x + GHOSTS; _k++) {
	point.i = _k;
#if dimension > 1
	for (point.j = GHOSTS; point.j < point.n.y + GHOSTS; point.j++)
#if dimension > 2
	  for (point.k = GHOSTS; point.k < point.n.z + GHOSTS; point.k++)
#endif
#endif
	    {...}
      }
  }
}

@define is_active(cell) (true)
@define is_leaf(cell)   (point.level == depth())
@define is_local(cell)  (true)
@define leaf            2
@def refine_cell(...) do {
  fprintf (stderr, "grid depths do not match. Aborting.\n");
  assert (0);
} while (0)
@
@define tree multigrid
#include "foreach_cell.h"

macro2 foreach_face_generic (char flags = 0, Reduce reductions = None,
			     const char * order = "xyz")
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    OMP(omp for schedule(static))
      for (_k = GHOSTS; _k <= point.n.x + GHOSTS; _k++) {
	point.i = _k;
#if dimension > 1
	for (point.j = GHOSTS; point.j <= point.n.y + GHOSTS; point.j++)
#if dimension > 2
	  for (point.k = GHOSTS; point.k <= point.n.z + GHOSTS; point.k++)
#endif
#endif
	    {...}
      }
  }
}

@define is_coarse() (point.level < depth())

#if dimension == 1
macro1 is_face_x() {
  {
    int ig = -1; NOT_UNUSED(ig);
    {...}
  }
}

// foreach_edge?

macro1 foreach_child (Point point = point, break = (_k = 2)) {
  {
    int _i = 2*point.i - GHOSTS;
    point.level++;
    point.n.x *= 2;
    for (int _k = 0; _k < 2; _k++) {
      point.i = _i + _k;
      POINT_VARIABLES();
      {...}
    }
    point.i = (_i + GHOSTS)/2;
    point.level--;
    point.n.x /= 2;
  }
}

#elif dimension == 2
#define foreach_edge() foreach_face(y,x)

macro1 is_face_x (Point p = point) {
  if (p.j < p.n.y + GHOSTS) {
    int ig = -1; NOT_UNUSED(ig);
    {...}
  }
}

macro1 is_face_y (Point p = point) {
  if (p.i < p.n.x + GHOSTS) {
    int jg = -1; NOT_UNUSED(jg);
    {...}
  }
}

macro1 foreach_child (Point point = point, break = (_k = _l = 2))
{
  {
    int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
    point.level++;
    point.n.x *= 2, point.n.y *= 2;
    for (int _k = 0; _k < 2; _k++)
      for (int _l = 0; _l < 2; _l++) {
	point.i = _i + _k; point.j = _j + _l;
	POINT_VARIABLES();
	{...}
      }
    point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2;
  }
}

#elif dimension == 3
macro1 is_face_x (Point p = point) {
  if (p.j < p.n.y + GHOSTS && p.k < p.n.z + GHOSTS) {
    int ig = -1; NOT_UNUSED(ig);
    {...}
  }
}

macro1 is_face_y (Point p = point) {
  if (p.i < p.n.x + GHOSTS && p.k < p.n.z + GHOSTS) {
    int jg = -1; NOT_UNUSED(jg);
    {...}
  }
}

macro1 is_face_z (Point p = point) {
  if (p.i < p.n.x + GHOSTS && p.j < p.n.y + GHOSTS) {
    int kg = -1; NOT_UNUSED(kg);
    {...}
  }
}

macro1 foreach_child (Point point = point, break = (_l = _m = _n = 2))
{
  {
    int _i = 2*point.i - GHOSTS;
    int _j = 2*point.j - GHOSTS;
    int _k = 2*point.k - GHOSTS;
    point.level++;
    point.n.x *= 2, point.n.y *= 2, point.n.z *= 2;
    for (int _l = 0; _l < 2; _l++)
      for (int _m = 0; _m < 2; _m++)
	for (int _n = 0; _n < 2; _n++) {
	  point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;
	  POINT_VARIABLES();
	  {...}
	}
    point.i = (_i + GHOSTS)/2;
    point.j = (_j + GHOSTS)/2;
    point.k = (_k + GHOSTS)/2;
    point.level--;
    point.n.x /= 2, point.n.y /= 2, point.n.z /= 2;
  }
}
#endif
  
@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

#include "neighbors.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  for (scalar s in list)
    if (!is_constant(s))
      for (int b = 0; b < s.block; b++) {
	real * data = (real *) multigrid->d;
	data += (s.i + b)*multigrid->shift[depth() + 1];
	for (size_t i = 0; i < multigrid->shift[depth() + 1]; i++, data++)
	  *data = val;
      }
}

// Boundaries

#if dimension == 1
macro2 foreach_boundary_dir (int l, int d, Reduce reductions = None)
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    SET_DIMENSIONS();
    if (d == left) {
      point.i = GHOSTS;
      ig = -1;
    }
    else if (d == right) {
      point.i = point.n.x + GHOSTS - 1;
      ig = 1;
    }
    {...}
  }
}

@define neighbor(o,p,q) ((Point){point.i+o, point.level, point.n _BLOCK_INDEX})
@define is_boundary(point) (point.i < GHOSTS || point.i >= point.n.x + GHOSTS)

#elif dimension == 2
macro2 foreach_boundary_dir (int l, int d, Reduce reductions = None)
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    SET_DIMENSIONS();
    int * _i = &point.j, _n = point.n.y;
    if (d == left) {
      point.i = GHOSTS;
      ig = -1;
    }
    else if (d == right) {
      point.i = point.n.x + GHOSTS - 1;
      ig = 1;
    }
    else if (d == bottom) {
      point.j = GHOSTS;
      _i = &point.i, _n = point.n.x;
      jg = -1;
    }
    else if (d == top) {
      point.j = point.n.y + GHOSTS - 1;
      _i = &point.i, _n = point.n.x;
      jg = 1;
    }
    int _l;
    OMP(omp for schedule(static))
      for (_l = 0; _l < _n + 2*GHOSTS; _l++) {
	*_i = _l;
	{...}
      }
  }
}

@def neighbor(o,p,q)
  ((Point){point.i+o, point.j+p, point.level, point.n _BLOCK_INDEX})
@
@def is_boundary(point) (point.i < GHOSTS || point.i >= point.n.x + GHOSTS ||
			 point.j < GHOSTS || point.j >= point.n.y + GHOSTS)
@

#elif dimension == 3
macro2 foreach_boundary_dir (int l, int d, Reduce reductions = None) {
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    SET_DIMENSIONS();
    int * _i = &point.j, * _j = &point.k;
    int _n[2] = { point.n.y, point.n.z };
    if (d == left) {
      point.i = GHOSTS;
      ig = -1;
    }
    else if (d == right) {
      point.i = point.n.x + GHOSTS - 1;
      ig = 1;
    }
    else if (d == bottom) {
      point.j = GHOSTS;
      _i = &point.i, _n[0] = point.n.x;
      jg = -1;
    }
    else if (d == top) {
      point.j = point.n.y + GHOSTS - 1;
      _i = &point.i, _n[0] = point.n.x;
      jg = 1;
    }
    else if (d == back) {
      point.k = GHOSTS;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = -1;
    }
    else if (d == front) {
      point.k = point.n.z + GHOSTS - 1;
      _i = &point.i; _j = &point.j;
      _n[0] = point.n.x, _n[1] = point.n.y;
      kg = 1;
    }
    int _l;
    OMP(omp for schedule(static))
      for (_l = 0; _l < _n[0] + 2*GHOSTS; _l++) {
	*_i = _l;
	for (int _m = 0; _m < _n[1] + 2*GHOSTS; _m++) {
	  *_j = _m;
	  {...}
	}
      }
  }
}

@def neighbor(o,p,q)
  ((Point){point.i+o, point.j+p, point.k+q, point.level, point.n _BLOCK_INDEX})
@
@def is_boundary(point) (point.i < GHOSTS || point.i >= point.n.x + GHOSTS ||
			 point.j < GHOSTS || point.j >= point.n.y + GHOSTS ||
			 point.k < GHOSTS || point.k >= point.n.z + GHOSTS)
@

#endif // dimension == 3

extern double (* default_scalar_bc[]) (Point, Point, scalar, bool *);
static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

macro2 foreach_boundary (int b, Reduce reductions = None)
{
  if (default_scalar_bc[b] != periodic_bc)
    foreach_boundary_dir (depth(), b, reductions)
      if (!is_boundary(point))
	{...}
}

@define neighborp(k,l,o) neighbor(k,l,o)

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int d = 0; d < 2*dimension; d++)
    if (default_scalar_bc[d] == periodic_bc)
      for (scalar s in scalars)
	if (!is_constant(s) && s.block > 0) {
	  if (is_vertex_scalar (s))
	    s.boundary[d] = s.boundary_homogeneous[d] = NULL;
	  else if (s.face) {
	    vector v = s.v;
	    v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
	  }
	}
  for (int bghost = 1; bghost <= BGHOSTS; bghost++)
    for (int d = 0; d < 2*dimension; d++) {

      scalar * list = NULL, * listb = NULL;
      for (scalar s in scalars)
	if (!is_constant(s) && s.block > 0) {
	  scalar sb = s;
#if dimension > 1
	  if (s.v.x.i >= 0) {
	    // vector component
	    int j = 0;
	    while ((&s.v.x)[j].i != s.i) j++;
	    sb = (&s.v.x)[(j - d/2 + dimension) % dimension];
	  }
#endif
	  if (sb.i >= 0 && sb.boundary[d] && sb.boundary[d] != periodic_bc) {
	    list = list_append (list, s);
	    listb = list_append (listb, sb);
	  }
	}
      
      if (list) {
	foreach_boundary_dir (l, d) {
	  scalar s, sb;
	  for (s,sb in list,listb) {
	    if ((s.face && sb.i == s.v.x.i) || is_vertex_scalar (s)) {
	      // normal component of face vector, or vertex scalar
	      if (bghost == 1)
		foreach_block()
		  s[(ig + 1)/2,(jg + 1)/2,(kg + 1)/2] =
		  sb.boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
	    }
	    else
	      // tangential component of face vector or centered
	      foreach_block()
		s[bghost*ig,bghost*jg,bghost*kg] =
		sb.boundary[d] (neighborp((1 - bghost)*ig,
					  (1 - bghost)*jg,
					  (1 - bghost)*kg),
				neighborp(bghost*ig,bghost*jg,bghost*kg),
				s, NULL);
	  }
	}
	free (list);
	free (listb);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

/* Periodic boundaries */

#if !_MPI

#if dimension == 1

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0 && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  if (l == 0) {
    foreach_level(0, noauto)
      for (scalar s in list1)
	for (int b = 0; b < s.block; b++) {
	  scalar sb = {s.i + b};
	  real v = sb[];
	  foreach_neighbor()
	    sb[] = v;
	}
    free (list1);
    return;
  }

  Point point = {0};
  point.level = l < 0 ? depth() : l;
  SET_DIMENSIONS();
  for (int i = 0; i < GHOSTS; i++)
    for (scalar s in list1)
      for (int b = 0; b < s.block; b++) {
	scalar sb = {s.i + b};
	sb[i] = sb[i + point.n.x];
      }
  for (int i = point.n.x + GHOSTS; i < point.n.x + 2*GHOSTS; i++)
    for (scalar s in list1)
      for (int b = 0; b < s.block; b++) {
	scalar sb = {s.i + b};
	sb[i] = sb[i - point.n.x];
      }
  free (list1);
}
    
#else // dimension != 1

@define VT _attribute[s.i].v.y

foreach_dimension()
static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0) {
      if (s.face) {
	scalar vt = VT;
	if (vt.i >= 0 && vt.boundary[right] == periodic_bc)
	  list1 = list_add (list1, s);
      }
      else if (s.boundary[right] == periodic_bc)
	list1 = list_add (list1, s);
    }
  if (!list1)
    return;

  if (l == 0) {
    foreach_level (0, noauto)
      for (scalar s in list1)
	for (int b = 0; b < s.block; b++) {
	  scalar sb = {s.i + b};
	  real v = sb[];
	  foreach_neighbor()
	    sb[] = v;
	}
    free (list1);
    return;
  }
  
  OMP_PARALLEL() {
    Point point = {0};
    point.level = l < 0 ? depth() : l;
    SET_DIMENSIONS();
#if dimension == 2  
    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n.y + 2*GHOSTS; j++) {
	for (int i = 0; i < GHOSTS; i++)
	  for (scalar s in list1)
	    for (int b = 0; b < s.block; b++) {
	      scalar sb = {s.i + b};
	      sb[i,j] = sb[i + point.n.x,j];
	    }
	for (int i = point.n.x + GHOSTS; i < point.n.x + 2*GHOSTS; i++)
	  for (scalar s in list1)
	    for (int b = 0; b < s.block; b++) {
	      scalar sb = {s.i + b};
	      sb[i,j] = sb[i - point.n.x,j];
	    }
      }
#else // dimension == 3
    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n.y + 2*GHOSTS; j++)
	for (int k = 0; k < point.n.z + 2*GHOSTS; k++) {
	  for (int i = 0; i < GHOSTS; i++)
	    for (scalar s in list1)
	      for (int b = 0; b < s.block; b++) {
		scalar sb = {s.i + b};
		sb[i,j,k] = sb[i + point.n.x,j,k];
	      }
	  for (int i = point.n.x + GHOSTS; i < point.n.x + 2*GHOSTS; i++)
	    for (scalar s in list1)
	      for (int b = 0; b < s.block; b++) {
		scalar sb = {s.i + b};
		sb[i,j,k] = sb[i - point.n.x,j,k];
	      }
	}
#endif
  }
  free (list1);
}

@undef VT

#endif // dimension != 1  
  
#endif // !_MPI

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = multigrid;
  free (m->d);
  free (m->shift);
  free (m);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}
 
void init_grid (int n)
{
  free_grid();
  Multigrid * m = qmalloc (1, Multigrid);
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n/Dimensions.x);
  N = (1 << grid->depth)*Dimensions.x;
#if !_MPI
  // mesh size
  grid->n = 1 << dimension*depth();
  foreach_dimension()
    grid->n *= Dimensions.x;
  grid->tn = grid->n;
#endif // !_MPI
  // box boundaries
  Boundary * b = qcalloc (1, Boundary);
  b->level = box_boundary_level;
  add_boundary (b);
#if _MPI
  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#else
  // periodic boundaries
  foreach_dimension() {
    Boundary * b = qcalloc (1, Boundary);
    b->level = periodic_boundary_level_x;
    add_boundary (b);
  }
#endif
  // allocate grid: this must be after mpi_boundary_new() since this modifies depth()
  m->shift = qmalloc (depth() + 2, size_t);
  size_t totalsize = 0;
  for (int l = 0; l <= depth() + 1; l++) {
    m->shift[l] = totalsize;
    struct _Point point = { .level = l };
    SET_DIMENSIONS();
    size_t size = 1;
    foreach_dimension()
      size *= point.n.x + 2*GHOSTS;
    totalsize += size;
  }
  m->d = (char *) malloc(m->shift[depth() + 1]*datasize);
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = multigrid;
  datasize += size;
  qrealloc (p->d, p->shift[depth() + 1]*datasize, char);
}

#if _MPI
int mpi_coords[dimension];
#undef _I
#undef _J
#undef _K
#define _I     (point.i - GHOSTS + mpi_coords[0]*(1 << point.level))
#define _J     (point.j - GHOSTS + mpi_coords[1]*(1 << point.level))
#define _K     (point.k - GHOSTS + mpi_coords[2]*(1 << point.level))
#endif

Point locate (double xp = 0, double yp = 0, double zp = 0)
{
  Point point = {0};
  point.level = depth();
  SET_DIMENSIONS();
  point.level = -1;
#if _MPI // fixme: can probably simplify with below
  point.i = (xp - X0)/L0*point.n.x*Dimensions.x + GHOSTS - mpi_coords[0]*point.n.x;
  if (point.i < GHOSTS || point.i >= point.n.x + GHOSTS)
    return point;
#if dimension >= 2
  point.j = (yp - Y0)/L0*point.n.x*Dimensions.x + GHOSTS - mpi_coords[1]*point.n.x;
  if (point.j < GHOSTS || point.j >= point.n.y + GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (zp - Z0)/L0*point.n.x*Dimensions.x + GHOSTS - mpi_coords[2]*point.n.x;
  if (point.k < GHOSTS || point.k >= point.n.z + GHOSTS)
    return point;
#endif  
#else // !_MPI
  point.i = (xp - X0)/L0*point.n.x + GHOSTS;
  if (point.i < GHOSTS || point.i >= point.n.x + GHOSTS)
    return point;
#if dimension >= 2
  point.j = (yp - Y0)/L0*point.n.x + GHOSTS;
  if (point.j < GHOSTS || point.j >= point.n.y + GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (zp - Z0)/L0*point.n.x + GHOSTS;
  if (point.k < GHOSTS || point.k >= point.n.z + GHOSTS)
    return point;
#endif  
#endif // !_MPI
  point.level = depth();
  return point;
}

#if _GPU
# include "variables.h"
#else
# include "multigrid-common.h"
#endif

macro2 foreach_vertex (char flags = 0, Reduce reductions = None) {
  foreach_face_generic (reductions) {
    int ig = -1; NOT_UNUSED(ig);
#if dimension > 1
    int jg = -1; NOT_UNUSED(jg);
#endif
#if dimension > 2
    int kg = -1; NOT_UNUSED(kg);
#endif
    {...}
  }
}

#if dimension == 3
macro foreach_edge (char flags = 0, Reduce reductions = None) {
  foreach_vertex (flags, reductions) {
    struct { int x, y, z; } _a = {point.i, point.j, point.k};
    foreach_dimension()
      if (_a.x < point.n.x + GHOSTS)
	{...}
  }
}
#endif // dimension == 3

ivec dimensions (int nx = 0, int ny = 0, int nz = 0)
{
  if (nx != 0 || ny != 0 || nz != 0) {
    Dimensions.x = Dimensions_scale = max(nx, 1);
#if dimension > 1
    Dimensions.y = max(ny, 1);
#endif
#if dimension > 2
    Dimensions.z = max(nz, 1);
#endif
  }
  return Dimensions;
}

#if _MPI
# include "multigrid-mpi.h"
#else // !_MPI

#define foreach_cell() foreach_cell_restore()

#endif // !_MPI
