typedef double number; // used in macros to indicate "any numeric type"
#define pi 3.14159265358979
#undef HUGE
#define HUGE 1e30f
#define nodata HUGE

macro number max (number a, number b) { return a > b ? a : b; }
macro number min (number a, number b) { return a < b ? a : b; }
macro number sq (number x) { return x*x; }
macro number cube (number x) { return x*x*x; }
macro int sign (number x) { return (int)(x > 0 ? 1 : -1); }
macro int sign2 (number x) { return (int)(x > 0 ? 1 : x < 0 ? -1 : 0); }
macro number clamp (number x, number a, number b) {
  return x < a ? a : x > b ? b : x;
}

#define swap(type,a,b) do { type _tmp_ = a; a = b; b = _tmp_; } while(false)

#include "grid/config.h"

static inline double noise() { return 1. - 2.*rand()/(double)RAND_MAX; }

// the grid
typedef struct {
  long n;       // number of (leaf) cells for this process
  long tn;      // number of (leaf) cells for all processes
  int depth;    // the depth for this process
  int maxdepth; // the maximum depth for all processes
} Grid;
Grid * grid = NULL;
// coordinates of the lower-left corner of the box
double X0 = 0., Y0 = 0., Z0 = 0.;
// size of the box
double L0 = 1. [1];
// number of grid points
#if dimension <= 2
int N = 64;
#else
int N = 16;
#endif

typedef struct { int i; } scalar;

typedef struct {
  scalar x;
#if dimension > 1
  scalar y;
#endif
#if dimension > 2
  scalar z;
#endif
} vector;

typedef struct {
  scalar * x;
#if dimension > 1
  scalar * y;
#endif
#if dimension > 2
  scalar * z;
#endif
} vectorl;

typedef struct {
  vector x;
#if dimension > 1
  vector y;
#endif
#if dimension > 2
  vector z;
#endif
} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
			   omp_out.x += omp_in.x,
			   omp_out.y += omp_in.y,
			   omp_out.z += omp_in.z))

#if dimension == 1
# define norm(v) fabs(v.x[])
# define dv() (Delta*cm[])
#elif dimension == 2
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[])))
# define dv() (sq(Delta)*cm[])
#else // dimension == 3
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[]) + sq(v.z[])))
# define dv() (cube(Delta)*cm[])
#endif

void normalize (coord * n)
{
  double norm = 0.;
  foreach_dimension()
    norm += sq(n->x);
  norm = sqrt(norm);
  foreach_dimension()
    n->x /= norm;
}

void origin (double x = 0., double y = 0., double z = 0.) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }

// boundary conditions for each direction/variable

#if dimension == 1
  enum { right, left };
#elif dimension == 2
  enum { right, left, top, bottom };
#else
  enum { right, left, top, bottom, front, back };
#endif
int nboundary = 2*dimension;

#define none -1

double  * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#include "grid/boundaries.h"

// attributes for each scalar

typedef struct {
  int x;
#if dimension > 1
  int y;
#endif
#if dimension > 2
  int z;
#endif
} ivec;
typedef double (* BoundaryFunc) (Point, Point, scalar, bool *);
typedef struct {
  BoundaryFunc * boundary;
  BoundaryFunc * boundary_homogeneous;
  double (* gradient)              (double, double, double);
  void   (* delete)                (scalar);
  char * name;
  ivec d; // staggering
  vector v;
  int face;
  bool   nodump, freed;
  int    block;
  scalar * depends; // boundary conditions depend on other fields
} _Attributes;

static _Attributes * _attribute = NULL;

#define foreach_block() // treatment of block data is disabled by default
#define foreach_blockf(s)

#if dimension == 1
ivec Dimensions = {1};
#elif dimension == 2
ivec Dimensions = {1,1};
#elif dimension == 3
ivec Dimensions = {1,1,1};
#endif

// lists

int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  qrealloc (list, len + 2, scalar);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  qrealloc (list, len + 2, scalar);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    for (scalar s1 in l)
      if (s1.i == s.i)
	return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    for (scalar s in l)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  for (scalar s in l2)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  for (scalar s in l)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", s.name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  qrealloc (list, len + 2, vector);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  for (vector w in list) {
    bool id = true;
    foreach_dimension()
      if (w.x.i != v.x.i)
	id = false;
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    for (vector v in l)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    foreach_dimension() {
      assert (s->i >= 0);
      v.x = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  for (tensor t in list) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  qrealloc (list, len + 2, tensor);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    foreach_dimension() {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  foreach_dimension()
    if (s.d.x != -1)
      return false;
  return true;
}

scalar * all = NULL; // all the fields
scalar * baseblock = NULL; // base block fields

// basic methods

scalar (* init_scalar)        (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector)        (vector, const char *);
vector (* init_face_vector)   (vector, const char *);
tensor (* init_tensor)        (tensor, const char *);
void   (* scalar_clone)       (scalar, scalar);

#define vector(x) (*((vector *)&(x)))

// timers

#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) + 
	  (tvend.tv_usec - t.tv.tv_usec)/1e6);
}

// Constant fields

const face vector zerof[] = {0.,0.,0.};
const face vector unityf[] = {1.,1.,1.};
const scalar unity[] = 1.;
const scalar zeroc[] = 0.;

// Metric

(const) face vector fm[] = {1.[0],1.[0],1.[0]};
(const) scalar cm[] = 1.[0];

// Embedded boundaries
// these macros are overloaded in embed.h

#define SEPS 0.
#define face_gradient_x(a,i) ((a[i] - a[i-1])/Delta)
#define face_gradient_y(a,i) ((a[0,i] - a[0,i-1])/Delta)
#define face_gradient_z(a,i) ((a[0,0,i] - a[0,0,i-1])/Delta)
#define face_value(a,i)      ((a[i] + a[i-1])/2.)
#define center_gradient(a)   ((a[1] - a[-1])/(2.*Delta))

// matrices

void * matrix_new (int n, int p, size_t size)
{
  void ** m = qmalloc (n, void *);
  char * a = qmalloc (n*p*size, char);
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (m[j][k]) >= big) {
	      big = fabs (m[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
	swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
	dum = m[ll][icol];
	m[ll][icol] = 0.0;
	for (l = 0; l < n; l++)
	  m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  free (((void **) m)[0]);
  free (m);
}

// Solver cleanup

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

// Default objects to display

static char * display_defaults = NULL;

static void free_display_defaults() {
  free (display_defaults);
}

void display (const char * commands, bool overwrite = false)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    free (display_defaults);
    display_defaults = malloc (strlen(commands) + 2);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = strdup ("@");
    display_defaults =
      realloc (display_defaults,
	       strlen(display_defaults) + strlen(commands) + 1);
    strcat (display_defaults, commands);
  }
}

#define display_control(val, ...)

typedef struct {
  double x;
#if dimension > 1
  double y;
#endif
#if dimension > 2
  double z;
#endif
} _coord;

// Types and macros for compatibility with GLSL

typedef struct {
  float r, g, b, a;
} vec4;

#if dimension == 1
# define avector(x, ...)    {x}
#elif dimension == 2
# define avector(x, y, ...) {x, y}
#else // dimension == 3
# define avector(x, y, z)   {x, y, z}
#endif

typedef struct {
  coord x, y, z;
} mat3;

OMP(omp declare reduction (+ : mat3 :
			   omp_out.x.x += omp_in.x.x,
			   omp_out.x.y += omp_in.x.y,
			   omp_out.x.z += omp_in.x.z,
			   omp_out.y.x += omp_in.y.x,
			   omp_out.y.y += omp_in.y.y,
			   omp_out.y.z += omp_in.y.z,
			   omp_out.z.x += omp_in.z.x,
			   omp_out.z.y += omp_in.z.y,
			   omp_out.z.z += omp_in.z.z
			   ))

typedef struct {
  uint32_t s;
} Adler32Hash;

static
inline void a32_hash_init (Adler32Hash * hash)
{
  hash->s = 0;
}

static
inline void a32_hash_add (Adler32Hash * hash, const void * data, size_t size)
{
  const uint8_t * buffer = (const uint8_t*) data;
  for (size_t n = 0; n < size; n++, buffer++)
    hash->s = *buffer + (hash->s << 6) + (hash->s << 16) - hash->s;
}

static
inline uint32_t a32_hash (const Adler32Hash * hash)
{
  return hash->s;
}
