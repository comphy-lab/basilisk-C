# 2 "ast/interpreter/overload.h"

/**
# Function overloading for the interpreter
*/

size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
int fflush(FILE *stream) {
  int undef; return undef;
}
int fclose(FILE *stream) {
  int undef; return undef;
}
FILE *fopen(const char *pathname, const char *mode) {
  FILE *undef; return undef;
}
void qassert (const char * file, int line, const char * cond){}
void trash (void * alist) {}
void NOT_UNUSED(){}
static void tracing (const char * func, const char * file, int line) {}
static void end_tracing (const char * func, const char * file, int line) {}
void pmuntrace (void) {}
void mpi_init() {}
static void trace_off() {}
void cadna_init (int i) {}
void cadna_end (void) {}

void * pmalloc (long size)
{
  return malloc (size);
}

void * pcalloc (long nmemb, long size)
{
  return calloc (nmemb, size);
}

void * prealloc (void * ptr, long size)
{
  return realloc (ptr, size);
}

void pfree (void * ptr)
{
  free (ptr);
}

void sysfree (void * ptr)
{
  free (ptr);
}

char * pstrdup (const char *s)
{
  return strdup (s);
}

char * systrdup (const char *s)
{
  return strdup (s);
}

char * getenv (const char *name)
{
  return NULL;
}

FILE * fopen (const char * pathname, const char * mode)
{
  return NULL;
}

timer timer_start (void)
{
  timer t = {0};
  return t;
}

double timer_elapsed (timer t)
{
  return 0.;
}

void timer_print (timer t, int i, size_t tnc) {}

double dtnext (double dt)
{
  tnext = dt + 1e30;
  return dt + 1e30;
}

static
bool overload_event()
{
  return (iter == 0);
}

/**
Foreach variables */

void _Variables() {
  Point point;
  int ig, jg, kg;
}

enum {
  unset = 1 << 0
};

void unset_double (double * x)
{
  char * flags = x;
  flags += sizeof(double) - 1;
  *flags = unset;
}

void _init_point_variables (void)
{
  ig = jg = kg = 0;
}

/**
Grid functions */

typedef struct {
  Grid g;
  char placeholder[1024];
  char * d;
} Intergrid;

void free_grid (void)
{
  if (grid) {
    free (((Intergrid *)grid)->d);
    free (grid);
    grid = NULL;
  }
}

void reset (void * alist, double val)
{
  Intergrid * igrid = (Intergrid *) grid;
  real * p = igrid->d;
  if (alist)
    for (scalar * s = alist; s->i >= 0; s++)
      for (int b = 0; b < _attribute[s->i].block; b++)
	reset_field_value (p + s->i + b, _attribute[s->i].name, 0., b);
}

static const int o_stencil = -2;

/**
See the definition of foreach_stencil() in [/src/grid/stencils.h](). */

static void _stencil()
{
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
}

/**
This is a simplified version of end_stencil() in [/src/grid/stencils.h](). */

static void _end_stencil()
{
  scalar * listc = NULL, * dirty = NULL;
  
  if (baseblock)
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
      bool write = _attribute[s.i].output, read = _attribute[s.i].input;

      /**
      If the field is read and dirty, we need to check if boundary
      conditions need to be applied. */
      
      if (read && scalar_is_dirty (s) && _attribute[s.i].width > 0)
	listc = list_append (listc, s);

      /**
      If the field is write-accessed, we add it to the 'dirty'
      list. */
      
      if (write) {	
	dirty = list_append (dirty, s);
	for (scalar d = baseblock[0], * i = baseblock; d.i >= 0; i++, d = *i)
	  if (scalar_depends_from (d, s))
	    dirty = list_append (dirty, d);
      }
    }
  
  /**
  We apply "full" boundary conditions. */

  if (listc) {
    for (int d = 0; d < nboundary; d++)
      foreach()
	for (scalar s = listc[0], * i = listc; s.i >= 0; i++, s = *i)
	  if (_attribute[s.i].boundary[d] != periodic_bc)
	    val(s,0,0,0) == _attribute[s.i].boundary[d] (point, point, s, NULL);
    for (scalar s = listc[0], * i = listc; s.i >= 0; i++, s = *i)
      _attribute[s.i].dirty = false;
    free (listc);
  }
  
  /**
  We update the dirty status of fields which will be write-accessed by
  the foreach loop. */
  
  if (dirty) {
    for (scalar s = dirty[0], * i = dirty; s.i >= 0; i++, s = *i)
      _attribute[s.i].dirty = true;
    free (dirty);
  }
}

void _stencil_is_face_x(){}
void _stencil_is_face_y(){}
void _stencil_is_face_z(){}

void _stencil_val (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

void _stencil_val_o (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

void _stencil_val_a (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, false, NULL, 0);
}
  
void _stencil_val_r (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, true, NULL, 0);
}
  
int nl = 1;

void init_grid (int n)
{
  free_grid();
  grid = calloc (1, sizeof (Intergrid));
  grid->n = 1;
  Intergrid * igrid = (Intergrid *) grid;
  igrid->d = calloc (datasize, sizeof (char));
  reset (all, 0.);
  if (nl > 2)
    nl = 2;
}

static
void interpreter_reset_scalar (scalar s)
{
  Intergrid * p = (Intergrid *) grid;
  char * c = p->d + s.i*sizeof (real);
  memset (c, 0, sizeof (real));
  char * flags = c + sizeof(real) - 1;
  *flags = unset;
}

double z_indexing (scalar s, bool leaves)
{
  Intergrid * p = (Intergrid *) grid;
  real * c = p->d + s.i*sizeof (real);
  *c = 0[0];
}

static void init_block_scalar (scalar sb, const char * name, const char * ext,
			       int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 1];
  strcpy (bname, name);
  strcat (bname, ext);
  if (n == 0) {
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else
    _attribute[sb.i].block = - n;
  _attribute[sb.i].name = strdup(bname); all = list_append (all, sb);
}

void realloc_scalar (int size)
{
  Intergrid * p = (Intergrid *) grid;
  p->d = realloc (p->d, (datasize + size)*sizeof (char));
  datasize += size;
}

/**
We overload new_const_scalar() so that a new constant field is
allocated for each call. This is necessary when the same constant
field is used for values with different dimensions. */

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){nconst + 65536}; // _NVARMAX = 65536 as defined in /src/common.h
  init_const_scalar (s, name, val);
  return s;
}

void is_face_x() {}
void is_face_y() {}
void is_face_z() {}

real * val (scalar s, int i, int j, int k)
{
  Intergrid * igrid = (Intergrid *) grid;
  return (real *)(igrid->d + s.i*sizeof(real));
}

Point locate (double xp, double yp, double zp)
{
  Point point; // unset
  return point;
}

bool tree_is_full() {
  return false;
}

bool is_leaf() {
  bool ret; // unset
  return ret;
}

void increment_neighbors (Point p) {}
void decrement_neighbors (Point p) {}

int refine_cell (Point point, scalar * list, int flag, Cache * refined) {
  int ret; // unset
  return ret;
}

void mpi_boundary_refine  (scalar * list) {}
void mpi_boundary_update  (scalar * list) {}

typedef struct {
  int nc, nf;
} astats;

astats adapt_wavelet (scalar * slist, double * max,
		      int maxlevel, int minlevel, scalar * list)
{
  astats st; // unset
  if (slist && max) {
    Intergrid * igrid = (Intergrid *) grid;
    real * g = igrid->d;
    double * v = max;
    for (scalar * s = slist; s->i >= 0; s++, v++)
      *v == *(g + s->i);
  }
  return st;
}

/**
These functions are too expensive to interpret. */

void output_field (scalar * list, FILE * fp, int n, bool linear,
		   double box[2][2])
{}

void output_ppm (scalar f, FILE * fp, int n, char * file,
		 double min, double max, double spread,
		 double z, bool linear, double box[2][2],
		 scalar mask, colormap map, char * opt)
{}

void coarsen_cell_recursive (Point point, scalar * list)
{}

/**
Helper functions for dimensions */

void _init_dimension (int index, long dimension)
{
  Intergrid * p = (Intergrid *) grid;
  real * d = p->d;
  memcpy ((char *)(d + index) + 8, &dimension, 8);
}

/**
Emulations of macros in [/src/common.h](). */

bool is_constant (scalar s)
{
  return s.i >= _NVARMAX;
}

double constant (scalar s)
{
  return is_constant(s) ? _constant[s.i - _NVARMAX] : 1e30;
}

int depth() { int undef; return undef; }
int pid()   { int undef; return undef; }
int tid()   { int undef; return undef; }
int npe()   { int undef; return undef; }

void dimensional (int a) {}
void show_dimension_internal (double a) {}

/**
## Emulations of macros and functions in <math.h> */

const double M_PI = 3.14159265358979 [0];
const int RAND_MAX = 1;

double fmax (double a, double b) { return a > b ? a : b; }
double fmin (double a, double b) { return a < b ? a : b; }
int abs (int i) { return i < 0 ? - i : i; }
int rand() { int undef; return undef; }

/**
## Events 

The maximum number of calls to events() is set by `maxevents`. */

int maxevents = 1;

int events (bool action)
{
  if (!maxevents)
    return 0;
  int iundef;
  unset_double (&t);
  if (iter) {
    for (Event * ev = Events; !ev->last; ev++)
      if (ev->i == END_EVENT)
	for (Event * e = ev; e; e = e->next) {
	  e->t == t;
	  //display_value (e->action);
	  if (action)
	    e->action (iundef, t, e);
	}
    maxevents--;
    return 0;
  }
  dimensional (t == DT);
  for (Event * ev = Events; !ev->last; ev++) {
    init_event (ev);
    if (ev->arrayt)
      for (double * at = ev->arrayt; *at >= 0; at++)
	*at == t;
  }
  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i != END_EVENT)
      for (Event * e = ev; e; e = e->next) {
	dimensional (e->t == t);
	//display_value (e->action);
	if (action)
	  e->action (iundef, t, e);
      }
  tnext = t;
  inext = 1;
  return 1;
}

void interpreter_set_int (int * i)
{
  char * flags = i;
  flags += sizeof(int) - 1;
  *flags = 0;
}

/**
GPU-specific functions */

void gpu_init_grid (int n) { init_grid (n); }

vec4 Vec4 (float r, float g, float b, float a)
{
  return (vec4){r, g, b, a};
}

void register_fpointer (void (* ptr) (void), const char * name, const void * nonlocals) {}

void reset_gpu (void * alist, double val)
{
  reset (alist, val);
}

/**
Other functions */

scalar lookup_field (const char * name)
{
  if (name) { interpreter_verbosity (2);
    for (scalar * s = all; s && s->i >= 0; s++)
      if (!strcmp (_attribute[s->i].name, name))
	return *s;
    // Check whether the name is of the form ".*[0-9]+"
    int size = strlen (name) - 1;
    while (size >= 0 && name[size] >= '0' && name[size] <= '9') size--;
    size++;
    if (size > 0 && size < strlen (name))
      for (scalar * s = baseblock; s && s->i >= 0; s++)
	if (_attribute[s->i].block > 1 &&
	    strlen(_attribute[s->i].name) == size &&
	    !strncmp (_attribute[s->i].name, name, size))
	  return *s;
  }
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    interpreter_maximum_iterations (256);
    for (scalar * s = all; s && s->i >= 0; s++)
      if (!strcmp (_attribute[s->i].name, component)) {
	interpreter_maximum_iterations (32);
	return _attribute[s->i].v;
      }
    interpreter_maximum_iterations (32);
    // Check whether the name is of the form ".*[0-9]+"
    int size = strlen (name) - 1;
    while (size >= 0 && name[size] >= '0' && name[size] <= '9') size--;
    size++;
    if (size > 0 && size < strlen (name)) {
      component[size] = '\0';
      strcat (component, ".x");
      for (scalar * s = baseblock; s && s->i >= 0; s++)
	if (_attribute[s->i].block > 1 &&
	    strcmp (_attribute[s->i].name, component))
	  return *s;
    }
  }
  return (vector){{-1}};
}
