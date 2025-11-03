#include "events.h"

void (* debug)    (Point);

@define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
@define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#include "fpe.h"
#include "stencils.h"

macro2 foreach_point (double _x = 0., double _y = 0., double _z = 0.,
		      char flags = 0, Reduce reductions = None)
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { _x, _y, _z };
    Point point = locate (_p.x, _p.y, _p.z); // fixme
    if (point.level >= 0)
      {...}
  }
}

macro2 foreach_region (coord p, coord box[2], coord n,
		       char flags = 0, Reduce reductions = None)
{
  {
    if (n.x < 1) n.x = 1;
    if (n.y < 1) n.y = 1;
    if (n.z < 1) n.z = 1;
    //  OMP(omp for schedule(static))
    for (int _i = 0; _i < (int) n.x; _i++) {
      p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);
      for (int _j = 0; _j < (int) n.y; _j++) {
	p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);
	for (int _k = 0; _k < (int) n.z; _k++) {
	  p.z = box[0].z + (box[1].z - box[0].z)/n.z*(_k + 0.5);
	  Point point = locate (p.x, p.y, p.z);
	  if (point.level >= 0) {
	    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
	    {...}
	  }
	}
      }
    }
  }
}

/**
Dirichlet and Neumann boundary conditions */

static inline
double dirichlet (double expr, Point point = point, scalar s = _s)
{
  return 2.*expr - s[];
}

static inline
double dirichlet_homogeneous (double expr, Point point = point, scalar s = _s)
{
  return - s[];
}

static inline
double dirichlet_face (double expr)
{
  return expr;
}

static inline
double dirichlet_face_homogeneous (double expr)
{
  return 0.;
}

static inline
double neumann (double expr, Point point = point, scalar s = _s)
{
  return Delta*expr + s[];
}

static inline
double neumann_homogeneous (double expr, Point point = point, scalar s = _s)
{
  return s[];
}

/**
Register functions on GPUs */

#if _GPU
#include "khash.h"
KHASH_MAP_INIT_INT64(PTR, External)

static khash_t(PTR) * _functions = NULL;
  
static External * _get_function (long ptr)
{
  if (!_functions)
    return NULL;  
  khiter_t k = kh_get (PTR, _functions, ptr);
  if (k == kh_end (_functions))
    return NULL;
  return &kh_value (_functions, k);
}

static void register_function (void (* ptr) (void), const char * name,
			       const char * kernel, const void * externals)
{
  static int index = 1;
  if (!_functions)
    _functions = kh_init (PTR), index = 1;
  int m = 0;
  for (const External * i = externals; i && i->name; i++, m++);
  External * copy = NULL;
  if (m > 0) {
    copy = malloc ((m + 1)*sizeof (External));
    memcpy (copy, externals, (m + 1)*sizeof (External));
  }
  int ret;
  khiter_t k = kh_put(PTR, _functions, (long) ptr, &ret);
  External p = {
    .name = (char *) name,
    .type = sym_function_definition,
    .pointer = (void *)(long) ptr, .nd = index++,
    .data = (void *) kernel, .externals = copy
  };
  kh_value(_functions, k) = p;
}

#define foreach_function(f, body) do {					\
    for (khiter_t k = kh_begin(_functions); k != kh_end(_functions); ++k) \
      if (kh_exist(_functions, k)) {					\
	External * f = &kh_value(_functions, k);			\
	body;								\
      }									\
  } while(0)

#endif // _GPU

/**
# Field allocation

If this routine is modified, do not forget to update [/src/ast/interpreter/overload.h](). */

static void init_block_scalar (scalar sb, const char * name, const char * ext,
			       int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    strcat (strcpy (bname, name), ext);
    sb.block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    sb.block = - n;
  }
  sb.name = strdup (bname);
  all = list_append (all, sb);
}

@define interpreter_set_int(...)
@define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(real);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && sb.freed)
      n++, sb.i++;
    if (n >= block) { // found n free slots
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
	init_block_scalar (sb, name, ext, n, block);
	interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}})); // fixme: only trashes one block?
      return s;
    }
    s.i = sb.i + 1;
  }
  
  // need to allocate new slots
  s = (scalar){nvar};
  assert (nvar + block <= _NVARMAX);

  if (_attribute == NULL)
    _attribute = (_Attributes *) calloc (nvar + block + 1, sizeof (_Attributes));
  else
    _attribute = (_Attributes *)
      realloc (_attribute, (nvar + block + 1)*sizeof (_Attributes));
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }
  // allocate extra space on the grid
  realloc_scalar (block*sizeof(real));
  trash (((scalar []){s, {-1}})); // fixme: only trashes one block?
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}
  
scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  foreach_dimension()
    v.x = alloc_block_scalar (name, ext.x, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    foreach_dimension()
      vb.x.i = v.x.i + i;
    init_vector (vb, NULL);
    foreach_dimension()
      vb.x.block = - i;
  }
  foreach_dimension()
    v.x.block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    foreach_dimension()
      vb.x.i = v.x.i + i;
    init_face_vector (vb, NULL);
    foreach_dimension()
      vb.x.block = - i;
  }
  foreach_dimension()
    v.x.block = block;  
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  tensor t;
  foreach_dimension() {
    strcat (strcpy (cname, name), ext.x);
    t.x = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  foreach_dimension()
    t.x.x = alloc_block_scalar (name, ext.x, 1);
  #if dimension > 1
    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;
  #endif
  #if dimension > 2
    t.x.z = alloc_block_scalar (name, ".x.z", 1);
    t.z.x = t.x.z;
    t.y.z = alloc_block_scalar (name, ".y.z", 1);
    t.z.y = t.y.z;
  #endif
  /* fixme: boundary conditions don't work!  This is because boundary
     attributes depend on the index and should (but cannot) be
     different for t.x.y and t.y.x */
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    qrealloc (_constant, s.i - _NVARMAX + 1, double);
    for (int i = nconst; i < s.i - _NVARMAX; i++)
      _constant[i] = 0.;
    nconst = s.i - _NVARMAX + 1;
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  foreach_dimension()
    init_const_scalar (v.x, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  foreach_dimension()
    v.x.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = clone.name;
  BoundaryFunc * boundary = clone.boundary;
  BoundaryFunc * boundary_homogeneous = clone.boundary_homogeneous;
  assert (src.block > 0 && clone.block == src.block);
  free (clone.depends);
  _attribute[clone.i] = _attribute[src.i];
  clone.name = cname;
  clone.boundary = boundary;
  clone.boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    clone.boundary[i] = src.boundary[i];
    clone.boundary_homogeneous[i] = src.boundary_homogeneous[i];
  }
  clone.depends = list_copy (src.depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(real), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  for (scalar s in l) {
    scalar c = s.block > 1 ? new_block_scalar("c", "", s.block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  for (scalar s in list)
    foreach_dimension()
      if (s.v.x.i >= 0 && map[s.v.x.i] >= 0)
	s.v.x.i = map[s.v.x.i];
  return list;
}

void delete (scalar * list)
{
  if (all == NULL) // everything has already been freed
    return;

  for (scalar f in list) {
    for (int i = 0; i < f.block; i++) {
      scalar fb = {f.i + i};
      if (f.delete)
	f.delete (fb);
      free (fb.name); fb.name = NULL;
      free (fb.boundary); fb.boundary = NULL;
      free (fb.boundary_homogeneous); fb.boundary_homogeneous = NULL;
      free (fb.depends); fb.depends = NULL;
      fb.freed = true;
    }
  }
  
  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  for (scalar f in list) {
    if (f.block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
	for (; s[f.block].i >= 0; s++)
	  s[0] = s[f.block];
	s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
	for (; s[1].i >= 0; s++)
	  s[0] = s[1];
	s->i = -1;
      }
    }
  }
}

void free_solver()
{  
  assert (_val_higher_dimension == 0.);

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }
  
  delete (all);
  free (all); all = NULL;
  free (baseblock); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      free (e);
      e = next;
    }
  }

  free (Events); Events = NULL;
  free (_attribute); _attribute = NULL;
  free (_constant); _constant = NULL;
#if _GPU
  foreach_function (f, free ((void *) f->externals));
  kh_destroy (PTR, _functions); _functions = NULL;
#endif
  free_grid();
  qpclose_all();
@if TRACE
  trace_off();
@endif
@if MTRACE
  pmuntrace();
@endif
@if _CADNA
  cadna_end();
@endif
}

// Cartesian methods

void (* boundary_level) (scalar *, int l);
void (* boundary_face)  (vectorl);

#define boundary(...)						\
  boundary_internal ((scalar *)__VA_ARGS__, __FILE__, LINENO)

void boundary_flux  (vector * list) __attribute__ ((deprecated));

void boundary_flux  (vector * list)
{
  vectorl list1 = {NULL};
  for (vector v in list)
    foreach_dimension()
      list1.x = list_append (list1.x, v.x);
  boundary_face (list1);
  foreach_dimension()
    free (list1.x);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t.i == s.i)
      return list;
  scalar * list1 = list;
  for (scalar d in _attribute[s.i].depends)
    if (d.dirty)
      list1 = list_add_depends (list1, d);
  return list_append (list1, s);
}

trace
void boundary_internal (scalar * list, const char * fname, int line)
{
  if (list == NULL)
    return;
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0) {
      if (scalar_is_dirty (s)) {
	if (s.face && s.dirty != 2)
	  foreach_dimension()
	    if (s.v.x.i == s.i)
	      listf.x = list_add (listf.x, s), flux = true;
	if (!is_constant(cm) && cm.dirty)
	  listc = list_add_depends (listc, cm);
	if (s.face != 2) // flux only
	  listc = list_add_depends (listc, s);
      }
#if 0
      else
	fprintf (stderr, "warning: bc already applied on '%s'\n", s.name);
#endif
    }
  if (flux) {
#if PRINTBOUNDARY
    int i = 0;
    foreach_dimension()
      if (listf.x) {
	fprintf (stderr, "boundary_internal: flux %c:", 'x' + i);
	for (scalar s in listf.x)
	  fprintf (stderr, " %d:%s", s.i, s.name);
	fputc ('\n', stderr);
      }
    i++;
#endif
    boundary_face (listf);
    foreach_dimension()
      free (listf.x);
  }
  if (listc) {
#if PRINTBOUNDARY
    fprintf (stderr, "boundary_internal: listc:");
    for (scalar s in listc)
      fprintf (stderr, " %d:%s", s.i, s.name);
    fputc ('\n', stderr);
#endif
    boundary_level (listc, -1);
    for (scalar s in listc)
      s.dirty = false;
    free (listc);
  }
}

void cartesian_boundary_level (scalar * list, int l)
{
  boundary_iterate (level, list, l);
}

void cartesian_boundary_face (vectorl vl)
{
  scalar * listc = NULL;
  foreach_dimension()
    for (scalar s in vl.x) {
      s.dirty = 2;
      listc = list_add_depends (listc, s);
    }
  boundary_level (listc, -1);
  free (listc);
}

static double symmetry (Point point, Point neighbor, scalar s, bool * data)
{
  return s[];
}

static double antisymmetry (Point point, Point neighbor, scalar s, bool * data)
{
  return -s[];
}

BoundaryFunc default_scalar_bc[] = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{
  // keep name
  char * pname;
  if (name) {
    free (s.name);
    pname = strdup (name);
  }
  else
    pname = s.name;
  int block = s.block;
  BoundaryFunc * boundary = s.boundary;
  BoundaryFunc * boundary_homogeneous = s.boundary_homogeneous;
  s.name = pname;
  if (block < 0)
    s.block = block;
  else
    s.block = block > 0 ? block : 1;
  /* set default boundary conditions */  
  s.boundary = boundary ? boundary : (BoundaryFunc *) malloc (nboundary*sizeof (BoundaryFunc));
  s.boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (BoundaryFunc *) malloc (nboundary*sizeof (BoundaryFunc));
  for (int b = 0; b < nboundary; b++)
    s.boundary[b] = s.boundary_homogeneous[b] =
      b < 2*dimension ? default_scalar_bc[b] : symmetry;
  s.gradient = NULL;
  foreach_dimension() {
    s.d.x = 0;  // not face
    s.v.x.i = -1; // not a vector component
  }
  s.face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  foreach_dimension()
    s.d.x = -1;
  for (int d = 0; d < nboundary; d++)
    s.boundary[d] = s.boundary_homogeneous[d] = NULL;
  return s;
}
  
BoundaryFunc default_vector_bc[] = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  foreach_dimension() {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    v.x.v = v;
  }
  /* set default boundary conditions */
  for (int d = 0; d < nboundary; d++)
    v.x.boundary[d] = v.x.boundary_homogeneous[d] =
      d < 2*dimension ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  foreach_dimension() {
    v.x.d.x = -1;
    v.x.face = true;
  }
  for (int d = 0; d < nboundary; d++)
    v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  foreach_dimension() {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  }
  /* set default boundary conditions */
  #if dimension == 1
    for (int b = 0; b < nboundary; b++)
      t.x.x.boundary[b] = t.x.x.boundary_homogeneous[b] =
	b < 2*dimension ? default_scalar_bc[b] : symmetry;
  #elif dimension == 2
    for (int b = 0; b < nboundary; b++) {
      t.x.x.boundary[b] = t.y.x.boundary[b] = 
	t.x.x.boundary_homogeneous[b] = t.y.y.boundary_homogeneous[b] = 
	b < 2*dimension ? default_scalar_bc[b] : symmetry;
      t.x.y.boundary[b] = t.y.y.boundary[b] = 
	t.x.y.boundary_homogeneous[b] = t.y.x.boundary_homogeneous[b] = 
	b < 2*dimension ? default_vector_bc[b] : antisymmetry;
    }
  #else
    assert (false); // not implemented yet
  #endif
  return t;
}

void output_cells (FILE * fp = stdout, coord c = {0}, double size = 0.)
{
  foreach() {
    bool inside = true;
    coord o = {x,y,z};
    foreach_dimension()
      if (inside && size > 0. &&
	  (o.x > c.x + size || o.x < c.x - size))
	inside = false;
    if (inside) {
      Delta /= 2.;
#if dimension == 1
      fprintf (fp, "%g 0\n%g 0\n\n", x - Delta, x + Delta);
#elif dimension == 2
      fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	       x - Delta, y - Delta,
	       x - Delta, y + Delta,
	       x + Delta, y + Delta,
	       x + Delta, y - Delta,
	       x - Delta, y - Delta);
#else // dimension == 3
      for (int i = -1; i <= 1; i += 2) {
	fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
		 x - Delta, y - Delta, z + i*Delta,
		 x - Delta, y + Delta, z + i*Delta,
		 x + Delta, y + Delta, z + i*Delta,
		 x + Delta, y - Delta, z + i*Delta,
		 x - Delta, y - Delta, z + i*Delta);
	for (int j = -1; j <= 1; j += 2)
	  fprintf (fp, "%g %g %g\n%g %g %g\n\n",
		   x + i*Delta, y + j*Delta, z - Delta,
		   x + i*Delta, y + j*Delta, z + Delta);
      }
#endif
    }
  }
  fflush (fp);
}

#if TREE && _MPI
static void output_cells_internal (FILE * fp)
{
  output_cells (fp);
}
#endif

static char * replace_ (const char * vname)
{
  char * name = strdup (vname), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
			const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp, 
	   "  load 'debug.plot'\n"
	   "  v=%s\n"
#if dimension == 1   
	   "  plot '%s' w l lc 0, "
	   "'%s' u 1+2*v:(0):2+2*v w labels tc lt 1 title columnhead(2+2*v)",
#elif dimension == 2
	   "  plot '%s' w l lc 0, "
	   "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",
#elif dimension == 3
	   "  splot '%s' w l lc 0, "
	   "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",
#endif
	   vname, cells, stencil);
  free (vname);
}

void cartesian_debug (Point point)
{
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  for (scalar v in all)
#if dimension == 1
    fprintf (fp, "x %s ", v.name);
#elif dimension == 2
    fprintf (fp, "x y %s ", v.name);
#elif dimension == 3
    fprintf (fp, "x y z %s ", v.name);
#endif
  fputc ('\n', fp);
  #if dimension == 1
    for (int k = -2; k <= 2; k++) {
      for (scalar v in all) {
	fprintf (fp, "%g ", x + k*Delta + v.d.x*Delta/2.);
	if (allocated(k))
	  fprintf (fp, "%g ", v[k]);
	else
	  fputs ("n/a ", fp);
      }
      fputc ('\n', fp);
    }
  #elif dimension == 2
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
	for (scalar v in all) {
	  fprintf (fp, "%g %g ",
		   x + k*Delta + v.d.x*Delta/2., 
		   y + l*Delta + v.d.y*Delta/2.);
	  if (allocated(k,l))
	    fprintf (fp, "%g ", v[k,l]);
	  else
	    fputs ("n/a ", fp);
	}
	fputc ('\n', fp);
      }
  #elif dimension == 3
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
	for (int m = -2; m <= 2; m++) {
	  for (scalar v in all) {
	    fprintf (fp, "%g %g %g ",
		     x + k*Delta + v.d.x*Delta/2., 
		     y + l*Delta + v.d.y*Delta/2.,
		     z + m*Delta + v.d.z*Delta/2.);
	    if (allocated(k,l,m))
	      fprintf (fp, "%g ", v[k,l,m]);
	    else
	      fputs ("n/a ", fp);
	  }
	  fputc ('\n', fp);
	}
  #endif
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp, 
	   "set term x11\n"
	   "set size ratio -1\n"
	   "set key outside\n");
  for (scalar s in all) {
    char * name = replace_ (s.name);
    fprintf (fp, "%s = %d\n", name, s.i);
    free (name);
  }
  fclose (fp);

  fprintf (stderr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (stderr, _attribute[0].name, name, stencil);
  fflush (stderr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar        = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector        = cartesian_init_vector;
  init_face_vector   = cartesian_init_face_vector;
  init_tensor        = cartesian_init_tensor;
  boundary_level     = cartesian_boundary_level;
  boundary_face      = cartesian_boundary_face;
  scalar_clone       = cartesian_scalar_clone;
  debug              = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
				  double xp = 0., double yp = 0., double zp = 0.)
{
#if dimension == 1
  x = (xp - x)/Delta - v.d.x/2.;
  int i = sign(x);
  x = fabs(x);
  /* linear interpolation */
  return v[]*(1. - x) + v[i]*x;
#elif dimension == 2
  x = (xp - x)/Delta - v.d.x/2.;
  y = (yp - y)/Delta - v.d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);
  /* bilinear interpolation */
  return ((v[]*(1. - x) + v[i]*x)*(1. - y) + 
	  (v[0,j]*(1. - x) + v[i,j]*x)*y);
#else // dimension == 3
  x = (xp - x)/Delta - v.d.x/2.;
  y = (yp - y)/Delta - v.d.y/2.;
  z = (zp - z)/Delta - v.d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);
  /* trilinear interpolation */
  return (((v[]*(1. - x) + v[i]*x)*(1. - y) + 
	   (v[0,j]*(1. - x) + v[i,j]*x)*y)*(1. - z) +
	  ((v[0,0,k]*(1. - x) + v[i,0,k]*x)*(1. - y) + 
	   (v[0,j,k]*(1. - x) + v[i,j,k]*x)*y)*z);
#endif  
}

trace
double interpolate (scalar v, double xp = 0., double yp = 0., double zp = 0.,
		    bool linear = true)
{
  double val = nodata;
  foreach_point (xp, yp, zp, reduction (min:val))
    val = linear ? interpolate_linear (point, v, xp, yp, zp) : v[];
  return val;
}

trace
void interpolate_array (scalar * list, coord * a, int n, double * v,
			bool linear = false)
{
  int len = 0;
  for (scalar s in list)
    len++;
  for (int i = 0; i < n; i++) {
    double * w = v;
#if _GPU
    coord p = a[i];
    for (scalar s in list) {
      double value = nodata;
      foreach_point (p.x, p.y, p.z, reduction(min:value))
	value = !linear ? s[] : interpolate_linear (point, s, p.x, p.y, 0.);
      *(w++) = value;
    }
#else // !_GPU
    for (scalar s in list)
      *(w++) = nodata;
    foreach_point (a[i].x, a[i].y, a[i].z, reduction(min:v[:len])) {
      int j = 0;
      for (scalar s in list)
	v[j++] = !linear ? s[] : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);
    }
#endif // !_GPU
    v = w;
  }
}

// Boundaries

typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  for (scalar s in all) {
    s.boundary = (BoundaryFunc *) realloc (s.boundary, nboundary*sizeof (BoundaryFunc));
    s.boundary_homogeneous =  (BoundaryFunc *)
      realloc (s.boundary_homogeneous, nboundary*sizeof (BoundaryFunc));
  }
  for (scalar s in all) {
    if (s.v.x.i < 0) // scalar
      s.boundary[b] = s.boundary_homogeneous[b] = symmetry;
    else if (s.v.x.i == s.i) { // vector
      vector v = s.v;
      foreach_dimension()
	v.y.boundary[b] = v.y.boundary_homogeneous[b] = symmetry;
      v.x.boundary[b] = v.x.boundary_homogeneous[b] =
	v.x.face ? NULL : antisymmetry;
    }
  }
  return b;
}

// Periodic boundary conditions

static double periodic_bc (Point point, Point neighbor, scalar s, bool * data)
{
  return s[];
}

static void periodic_boundary (int d)
{
  /* We change the conditions for existing scalars. */
  for (scalar s in all)
    if (is_vertex_scalar (s))
      s.boundary[d] = s.boundary_homogeneous[d] = NULL;
    else
      s.boundary[d] = s.boundary_homogeneous[d] = periodic_bc;
  /* Normal components of face vector fields should remain NULL. */
  for (scalar s in all)
    if (s.face) {
      vector v = s.v;
      v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
    }
  /* We also change the default boundary conditions (for new fields). */
  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{
  #if dimension < 2
    assert (dir <= left);
  #elif dimension < 3
    assert (dir <= bottom);
  #else
    assert (dir <= back);
  #endif
  // This is the component in the given direction i.e. 0 for x and 1 for y
  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}

// for debugging
double getvalue (Point point, scalar s, int i, int j, int k)
{
  return s[i,j,k];
}

void default_stencil (Point p, scalar * list)
{
  for (scalar s in list) {
    if (s.v.x.i != -1) {
      vector v = s.v;
      for (scalar c in {v})
	c.input = c.output = c.nowarning = true, c.width = 2;
    }
    else
      s.input = s.output = s.nowarning = true, s.width = 2;
  }
}

/**
This displays a (1D,2D,3D) stencil index. */

static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < dimension; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
		  const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  if (s.block < 0)
    s.i += s.block;
  if (!s.name) {
    fprintf (stderr, "%s:%d: error: trying to access a deleted field\n",
	     file, line);
    exit (1);
  }
  int index[] = {i, j, k};
  for (int d = 0; d < dimension; d++) {
    if (index[d] == o_stencil)
      index[d] = 2;
    else
      index[d] += (&p.i)[d];
  }
  bool central = true;
  for (int d = 0; d < dimension; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
	       file, line, s.name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!s.output)
      s.input = true;
  }
  else {
    s.input = true;
    int d = 0;
    foreach_dimension() {
      if ((!s.face || s.v.x.i != s.i) && abs(index[d]) > s.width)
	s.width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
		    const char * file, int line)
{
  if (is_constant(s) || s.i < 0) {
    fprintf (stderr, "%s:%d: error: trying to modify a%s field\n",
	     file, line, s.i < 0 ? "n undefined" : " constant");
    exit (1);
  }
  if (s.block < 0)
    s.i += s.block;
  if (!s.name) {
    fprintf (stderr, "%s:%d: error: trying to access a deleted field\n",
	     file, line);
    exit (1);
  }
  int index[] = {i, j, k};
  for (int d = 0; d < dimension; d++)
    index[d] += (&p.i)[d];    
  for (int d = 0; d < dimension; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
	       file, line, s.name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !s.output)
    s.input = true;
  s.output = true;
}
