#define GRIDNAME "Cartesian"
#define dimension 2
#define GHOSTS 1

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

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

/*
  The OpenACC version currently uses a simple pointer "acc_dataarr"
*/
@def data(a,k,l,m) acc_dataarr[(point.n + 2)*(point.n + 2)*a.i +
			       (point.i + k)*(point.n + 2) +
			       (point.j + l)] @

/*
  These global variables are defined in accelerator memory
  and updated before every "foreach"-type loop runs.
 */
int N;
double X0, Y0, L0, G, CFL;
ACC(acc declare create(N, X0, Y0, L0, G, CFL))

@define allocated(...) true

@define POINT_VARIABLES VARIABLES

// --------------------------------------------------------------------------------

// Accelerator version for "foreach"
/*
  Note: OpenACC currently works best when simple variables types are used,
  the implementation thus uses "iacc" and "jacc" as loop variables, and
  copies "ngrid" across as a scalar variable. Variable "point" needs to
  be thread private and is initialised inside the loop.
*/
@def foreach_acc(clause)
  _OMPSTART // Needed for parallel reduction treatment in qcc
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  size_t ngrid = cartesian->n;
  size_t nall = list_len(all)+1; NOT_UNUSED(nall);
  ACC(acc update device(N, X0, Y0, L0, G, CFL))
  ACC(acc parallel loop gang clause independent present(acc_dataarr) copyin(all[0:nall])) // NEED GANG AND VECTOR HERE???
  for (size_t iacc = 1; iacc <= ngrid; iacc++) {
    ACC(acc loop vector independent)
    for (size_t jacc = 1; jacc <= ngrid; jacc++) {
      Point point;
      point.n = ngrid;
      point.i = iacc;
      point.j = jacc;
      POINT_VARIABLES
@
@define end_foreach_acc() }} _OMPEND

// Accelerator version for "foreach_face"
@def foreach_face_acc_generic(clause)
  _OMPSTART // Needed for parallel reduction treatment in qcc
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  size_t ngrid = cartesian->n;
  size_t nall = list_len(all)+1; NOT_UNUSED(nall);
  ACC(acc update device(N, X0, Y0, L0, G, CFL))
  ACC(acc parallel loop gang clause independent present(acc_dataarr) copyin(all[0:nall]))
  for (size_t iacc = 1; iacc <= ngrid + 1; iacc++) {
    ACC(acc loop vector independent)
    for (size_t jacc = 1; jacc <= ngrid + 1; jacc++) {
      Point point;
      point.n = ngrid;
      point.i = iacc;
      point.j = jacc;
      POINT_VARIABLES
@
@define end_foreach_face_acc_generic() }} _OMPEND

@def foreach_vertex_acc()
foreach_face_acc_generic() {
  x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex_acc() } end_foreach_face_acc_generic()

// --------------------------------------------------------------------------------

@def foreach(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 1; _k <= point.n; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n; point.j++) {
      POINT_VARIABLES
@
@define end_foreach() }} OMP_END_PARALLEL()

@def foreach_face_generic(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 1; _k <= point.n + 1; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n + 1; point.j++) {
      POINT_VARIABLES
@
@define end_foreach_face_generic() }} OMP_END_PARALLEL()

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

#define foreach_edge() foreach_face(y,x)

@define is_face_x() (point.j <= point.n)
@define is_face_y() (point.i <= point.n)

@if TRASH
@ undef trash
@ define trash cartesian_trash
@endif

#include "neighbors.h"

void cartesian_trash (void * alist)
{
  scalar * list = alist;
  size_t np = acc_dataarr_len/datasize*sizeof(double);
  for (scalar s in list) {
    if (!is_constant(s)) {
      ACC(acc kernels loop present(acc_dataarr))
      for (size_t i = 0; i < np; i++) {
	acc_dataarr[s.i*np + i] = undefined;
      }
    }
  }
}

// ghost cell coordinates for each direction
static int _ig[] = {1,-1,0,0}, _jg[] = {0,0,1,-1};

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  Point point;
  point.n = cartesian->n;

  if (d % 2)
    ig = jg = 0;
  else {
    ig = _ig[d]; jg = _jg[d];
  }
  int _start = GHOSTS, _end = point.n + GHOSTS;  

  for (scalar s in list) {
    scalar b = s.v.x;
    b.boundary[d] (s, _start, _end, d, ig, jg);
  }
}

static void box_boundary_level_tangent (const Boundary * b, 
					scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  Point point;
  point.n = cartesian->n;

  ig = _ig[d]; jg = _jg[d];
  int _start = GHOSTS, _end = point.n + 2*GHOSTS;

  for (scalar s in list) {
    scalar b = s.v.y;
    b.boundary[d] (s, _start, _end, d, ig, jg);
  }
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL, * tangent = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
	else {
	  scalar b = s.v.y;
	  if (b.boundary[d])
	    tangent = list_add (tangent, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  Point point;
  point.n = cartesian->n;

  /* 
     Boundary enumeration:
     d=0  right boundary
     d=1  left boundary
     d=2  top boundary
     d=3  bottom boundary

     Index increments for ghost zones:
     _ig = {1, -1, 0,  0} for right/left/top/bottom
     _jg = {0,  0, 1, -1} for right/left/top/bottm
  */

  ig = _ig[d]; jg = _jg[d]; // Get ghost zone index steps (0 <= d < 4 => denotes boundary)
  int _start = 1, _end = point.n;  // Start and stop for boundary loop, includes corners
  /* traverse corners only for top and bottom */
  if (d > left) { _start--; _end++; } // Do not compute corners for top and bottom boundaries

  for (scalar s in centered) { // loop over all centered fields
    scalar b = (s.v.x.i < 0 ? s :
		s.i == s.v.x.i && d < top ? s.v.x :
		s.i == s.v.y.i && d >= top ? s.v.x :
		s.v.y);
    b.boundary[d] (s, _start, _end, d, ig, jg);
  }

  free (centered);

  box_boundary_level_normal (b, normal, l);
  free (normal);
  box_boundary_level_tangent (b, tangent, l);
  free (tangent);

}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();

  /* Free acc_dataarr on accelerator */
  ACC(acc exit data delete(acc_dataarr[0:acc_dataarr_len]))
  
  /* Point acc_datarr to NULL on host */
  acc_dataarr = NULL;

  /* Set size of array to zero */
  acc_dataarr_len = 0;

  free (cartesian->d);
  free (cartesian);
  grid = NULL;
}

void init_grid (int n)
{
  if (cartesian && n == cartesian->n)
    return;
  free_grid();
  Cartesian * p = malloc(sizeof(Cartesian));
  size_t len = (n + 2)*(n + 2)*datasize;
  p->n = N = n;
  p->d = malloc (len);

  /* Point acc_dataarr to host data array and store size */
  acc_dataarr = (double *)p->d;
  acc_dataarr_len = len/sizeof(double);

  /* Create acc_dataarr on accelerator and set
     values to undefined.  */
  ACC(acc enter data create(acc_dataarr[0:acc_dataarr_len]))
  ACC(acc kernels loop present(acc_dataarr))
  for (size_t i = 0; i < len/sizeof(double); i++)
    acc_dataarr[i] = undefined;

  grid = (Grid *) p;
  trash (all);
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
    add_boundary (b);
  }
  // mesh size
  grid->n = grid->tn = sq(n);
}

void realloc_scalar (void)
{
  Cartesian * p = cartesian;
  size_t len = (p->n + 2)*(p->n + 2);

  /* Copy to host and free array on accelerator */
  ACC(acc exit data copyout(acc_dataarr[0:acc_dataarr_len]))

  /* Reset array pointer */
  acc_dataarr = NULL;

  p->d = realloc (p->d, len*datasize);

  /* Set array pointer and size */
  acc_dataarr = (double *)p->d;
  acc_dataarr_len = len*datasize/sizeof(double);

  /* Create data array on accelerator and copy values */
  ACC(acc enter data copyin(acc_dataarr[0:acc_dataarr_len]))

}

/*
  Accelerator version of "locate" - returns "point" in
  function argument rather than as return value, due to
  the latter leading to compiler failures with PGI 16.1.
*/
ACC(acc routine seq)
void locate2Dacc (double x, double y, Point * point)
{
  point -> n = N;
  point -> i = (x - X0)/L0*(point -> n) + 1;
  point -> j = (y - Y0)/L0*(point -> n) + 1;
  point -> level = (point -> i >= 1 && point -> i <= point -> n &&
		    point -> j >= 1 && point -> j <= point -> n) ? 0 : - 1;
}

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point;
  point.n = cartesian->n;
  point.i = (p.x - X0)/L0*point.n + 1;
  point.j = (p.y - Y0)/L0*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n &&
		 point.j >= 1 && point.j <= point.n) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"

// Call "cartesian_methods" to avoid code replication
void cartesianacc_methods() {
  cartesian_methods();
}
