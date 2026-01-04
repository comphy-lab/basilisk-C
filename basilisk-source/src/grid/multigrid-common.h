#define MULTIGRID 1

#include "variables.h"
#include "cartesian-common.h"

auto macro2 foreach_level_or_leaf (int l, char flags = 0, Reduce reductions = None)
{
  foreach_level (l, flags, reductions)
    {...}
}

auto macro2 foreach_coarse_level (int l, char flags = 0, Reduce reductions = None)
{
  foreach_level (l, flags, reductions)
    {...}
}

// scalar attributes

attribute {
  void (* prolongation) (Point, scalar);
  void (* restriction)  (Point, scalar);
}

// Multigrid methods

void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{
  foreach_blockf (s) {
    double sum = 0.;
    foreach_child()
      sum += s[];
    s[] = sum/(1 << dimension);
  }
}

static inline void restriction_volume_average (Point point, scalar s)
{
  foreach_blockf (s) {
    double sum = 0.;
    foreach_child()
      sum += cm[]*s[];
    s[] = sum/(1 << dimension)/(cm[] + 1e-30);
  }
}

static inline void face_average (Point point, vector v)
{
  foreach_dimension()
    foreach_blockf (v.x) {
    #if dimension == 1
      v.x[] = fine(v.x,0);
      v.x[1] = fine(v.x,2);
    #elif dimension == 2
      v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
      v.x[1] = (fine(v.x,2,0) + fine(v.x,2,1))/2.;
    #else // dimension == 3
      v.x[] = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
	       fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      v.x[1] = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
		fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;
    #endif
  }
}
  
static inline void restriction_face (Point point, scalar s)
{
  face_average (point, s.v);
}

static inline void restriction_vertex (Point point, scalar s)
{
  foreach_blockf (s)
    for (int i = 0; i <= 1; i++) {
      s[i] = fine(s,2*i);
#if dimension >= 2  
      s[i,1] = fine(s,2*i,2);
#endif
#if dimension >= 3
      for (int j = 0; j <= 1; j++)
	s[i,j,1] = fine(s,2*i,2*j,2);
#endif
    }
}

static inline void no_restriction (Point point, scalar s) {}

static inline void no_data (Point point, scalar s) {
  foreach_child()
    foreach_blockf (s)
      s[] = nodata;
}

void wavelet (scalar s, scalar w)
{
  restriction ({s});
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    foreach_coarse_level (l, nowarning) {
      foreach_child()
        w[] = s[];
      s.prolongation (point, s);
      foreach_child() {
        double sp = s[];
        s[] = w[];
        /* difference between fine value and its prolongation */
        w[] -= sp;
      }
    }
    boundary_level ({w}, l + 1);
  }
  /* root cell */
  foreach_level(0) 
    w[] = s[];
  boundary_level ({w}, 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  foreach_level(0) 
    s[] = w[];
  boundary_level ({s}, 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    foreach_coarse_level (l, nowarning) {
      s.prolongation (point, s);
      foreach_child()
        s[] += w[];
    }
    boundary_level ({s}, l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{
  #if dimension == 1
    return (3.*coarse(s) + coarse(s,child.x))/4.;
  #elif dimension == 2
    return (9.*coarse(s) + 
	    3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
	    coarse(s,child.x,child.y))/16.;
  #else // dimension == 3
    return (27.*coarse(s) + 
	    9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		coarse(s,0,0,child.z)) + 
	    3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		coarse(s,0,child.y,child.z)) + 
	    coarse(s,child.x,child.y,child.z))/64.;
  #endif
}

static inline void refine_bilinear (Point point, scalar s)
{
  foreach_child()
    foreach_blockf (s)
      s[] = bilinear (point, s);
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{
#if dimension == 1
  return quadratic (coarse(s,0), coarse(s,child.x), coarse(s,-child.x));
#elif dimension == 2
  return
    quadratic (quadratic (coarse(s,0,0),
			  coarse(s,child.x,0),
			  coarse(s,-child.x,0)),
	       quadratic (coarse(s,0,child.y),
			  coarse(s,child.x,child.y),
			  coarse(s,-child.x,child.y)),
	       quadratic (coarse(s,0,-child.y),
			  coarse(s,child.x,-child.y),
			  coarse(s,-child.x,-child.y)));
#else // dimension == 3
  assert (false);
  return 0.;
#endif
}

static inline double biquadratic_vertex (Point point, scalar s)
{
#if dimension == 1
  return (6.*s[] + 3.*s[-1] - s[1])/8.;
#elif dimension == 2
  return (36.*s[] + 18.*(s[-1] + s[0,-1]) - 6.*(s[1] + s[0,1]) +
	  9.*s[-1,-1] - 3.*(s[1,-1] + s[-1,1]) + s[1,1])/64.;  
#elif dimension == 3
  assert (false);
  return 0.;  
#endif
}

static inline void refine_biquadratic (Point point, scalar s)
{
  foreach_child()
    foreach_blockf (s)
      s[] = biquadratic (point, s);
}

static inline void refine_linear_single (Point point, scalar s)
{
  coord g;
  if (s.gradient)
    foreach_dimension()
      g.x = s.gradient (s[-1], s[], s[1]);
  else
    foreach_dimension()
      g.x = (s[1] - s[-1])/2.;

  double sc = s[], cmc = 4.*cm[], sum = cm[]*(1 << dimension);
  foreach_child() {
    s[] = sc;
    foreach_dimension()
      s[] += child.x*g.x*cm[-child.x]/cmc;
    sum -= cm[];
  }
  assert (fabs(sum) < 1e-10);
}

static inline void refine_linear (Point point, scalar s)
{
  foreach_blockf (s)
    refine_linear_single (point, s);
}

static inline void refine_reset (Point point, scalar v)
{
  foreach_child()
    foreach_blockf (v)
      v[] = 0.;
}

static inline void refine_injection (Point point, scalar v)
{
  double val = v[];
  foreach_child()
    foreach_blockf (v) 
      v[] = val;
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  s.prolongation = refine_bilinear;
  s.restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  s.restriction = restriction_vertex;
  return s;
}

static void multigrid_setup_vector (vector v)
{
  foreach_dimension() {
    v.x.prolongation = refine_bilinear;
    v.x.restriction = restriction_average;
  }  
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  foreach_dimension()
    v.y.restriction = no_restriction;
  v.x.restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  foreach_dimension()
    multigrid_setup_vector (t.x);
  return t;
}

void multigrid_debug (Point point)
{
  cartesian_debug (point);
  
  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
    #if dimension == 1
      double xc = x - child.x*Delta/2.;
      for (int k = 0; k <= 1; k++)
	for (scalar v in all)
	  fprintf (fp, "%g %g ", 
		   xc + k*child.x*Delta*2. + v.d.x*Delta, 
		   coarse(v,k*child.x));
      fputc ('\n', fp);
      fprintf (stderr, ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 3 t ''", name);
      fprintf (plot,   ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 3 t ''", name);
    #elif dimension == 2
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
	for (int l = 0; l <= 1; l++) {
	  for (scalar v in all)
	    fprintf (fp, "%g %g %g ", 
		     xc + k*child.x*Delta*2. + v.d.x*Delta, 
		     yc + l*child.y*Delta*2. + v.d.y*Delta,
		     coarse(v,k*child.x,l*child.y));
	  fputc ('\n', fp);
	}
      fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot,   ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
    #elif dimension == 3
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
	for (int l = 0; l <= 1; l++)
	  for (int m = 0; m <= 1; m++) {
	    for (scalar v in all)
	      fprintf (fp, "%g %g %g %g ", 
		       xc + k*child.x*Delta*2. + v.d.x*Delta, 
		       yc + l*child.y*Delta*2. + v.d.y*Delta,
		       zc + m*child.z*Delta*2. + v.d.z*Delta,
		       coarse(v,k*child.x,l*child.y,m*child.z));
	    fputc ('\n', fp);
	  }
      fprintf (stderr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
	       name);
      fprintf (plot,   ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
	       name);
    #endif
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
    #if dimension == 1
      double xf = x - Delta/4.;
      for (int k = -2; k <= 3; k++)
	for (scalar v in all) {
	  fprintf (fp, "%g ", xf + k*Delta/2. + v.d.x*Delta/4.);
	  if (allocated_child(k))
	    fprintf (fp, "%g ", fine(v,k));
	  else
	    fputs ("n/a ", fp);
	}
      fputc ('\n', fp);
      fprintf (stderr, ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 2 t ''", name);
      fprintf (plot,   ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 2 t ''", name);
    #elif dimension == 2
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
	for (int l = -2; l <= 3; l++) {
	  for (scalar v in all) {
	    fprintf (fp, "%g %g ", 
		     xf + k*Delta/2. + v.d.x*Delta/4., 
		     yf + l*Delta/2. + v.d.y*Delta/4.);
	    if (allocated_child(k,l))
	      fprintf (fp, "%g ", fine(v,k,l));
	    else
	      fputs ("n/a ", fp);
	  }
	  fputc ('\n', fp);
	}
      fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot,   ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
    #elif dimension == 3
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
	for (int l = -2; l <= 3; l++)
	  for (int m = -2; m <= 3; m++) {
	    for (scalar v in all) {
	      fprintf (fp, "%g %g %g ", 
		       xf + k*Delta/2. + v.d.x*Delta/4., 
		       yf + l*Delta/2. + v.d.y*Delta/4.,
		       zf + m*Delta/2. + v.d.z*Delta/4.);
	      if (allocated_child(k,l,m))
		fprintf (fp, "%g ", fine(v,k,l,m));
	      else
		fputs ("n/a ", fp);
	    }
	    fputc ('\n', fp);
	  }
      fprintf (stderr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
	       name);
      fprintf (plot,   ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
	       name);
    #endif
    fclose (fp);
  }
  fflush (stderr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  for (scalar s in list) 
    if (!is_constant (s) && s.block > 0) {
      if (s.restriction == restriction_average) {
	listdef = list_add (listdef, s);
	list2 = list_add (list2, s);
      }
      else if (s.restriction != no_restriction) {
	listc = list_add (listc, s);
	if (s.face)
	  foreach_dimension()
	    list2 = list_add (list2, s.v.x);
	else
	  list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      foreach_coarse_level(l, nowarning) {
	for (scalar s in listdef)
	  restriction_average (point, s);
	for (scalar s in listc) {
	  s.restriction (point, s);
	}
      }
      boundary_iterate (level, list2, l);      
    }
    free (listdef);
    free (listc);
    free (list2);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar        = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector        = multigrid_init_vector;
  init_face_vector   = multigrid_init_face_vector;
  init_tensor        = multigrid_init_tensor;
  restriction        = multigrid_restriction;
  debug              = multigrid_debug;
}

/**
## Size of subtrees

The function below store in *size* the number of cells (or leaves if
*leaves* is set to *true*) of each subtree. */

void subtree_size (scalar size, bool leaves)
{

  /**
  The size of leaf "subtrees" is one. */

  foreach()
    size[] = 1;
  
  /**
  We do a (parallel) restriction to compute the size of non-leaf
  subtrees. */

  boundary_iterate (restriction, {size}, depth());
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level(l) {
      double sum = !leaves;
      foreach_child()
	sum += size[];
      size[] = sum;
    }
    boundary_iterate (restriction, {size}, l);
  }
}
