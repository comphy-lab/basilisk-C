/**
# A range of GPU tests */

#include "grid/multigrid.h"
#include "utils.h"

scalar s[], s1[];
vector v[];

void variable_list (scalar * list)
{
  foreach() {
    int i = 0;
    for (scalar s in list)
      s[] = 13 + i++;
  }
}

double myfunc2 (double x)
{
  return x;
}

double myfunc (double x)
{
  return myfunc2 (x);
}

void myfunc3 (double * a)
{
  *a += 2;
}

void myfunc4 (Point point, scalar a, (const) scalar b)
{
  a[] = b[];
}

double myfunc5 (double x)
{
  return x*x;
}

void myfunc6 (Point point, scalar s)
{
  s[] = v.x[];
}

double vv = 3;

void myfunc7 (Point point, scalar s, double a)
{
  s[] = a/vv;
}

void myfunc8 (scalar s)
{
  foreach()
    s[] = 0;
}

attribute {
  double (* func) (double x);
}

int main (int argc, char * argv[])
{
  size (1[0]);
  init_grid (1);

#if 0  
  gpu_limits (stderr);
#endif

  /**
  ## reset() */

  reset ({s}, 1);
  reset ({s1}, 2);
  foreach (serial)
    fprintf (stderr, "0) %g %g\n", s[], s1[]);
  
  /**
  ## Check input/output for vector fields */

  foreach()
    v.x[] = 1., v.y[] = 2.;
  foreach (serial)
    fprintf (stderr, "1) v.x: %g v.y: %g\n", v.x[], v.y[]);
  foreach()
    v.x[] += 1., v.y[] += 2.;
  foreach (serial)
    fprintf (stderr, "2) v.x: %g v.y: %g\n", v.x[], v.y[]);

  /**
  ## Check consistent writes to individual texture components */
  
  foreach()
    v.y[] = 3; // v.x[] should not be modified
  foreach (serial)
    fprintf (stderr, "3) v.x: %g v.y: %g\n", v.x[], v.y[]);

  /**
  ## Check consistent CPU copies of individual texture components */
  
  foreach (cpu)
    v.x[] = 4.;
  foreach()
    v.x[] += v.y[];
  foreach (serial)
    fprintf (stderr, "4) v.x: %g v.y: %g\n", v.x[], v.y[]);

  /**
  ## Check for "no inputs" */

  {
    scalar a[], b[];
    foreach() {
      a[] = 1.;
      b[] = 3. - a[];
    }
    foreach (serial)
      fprintf (stderr, "5) %g %g\n", a[], b[]);
  }
  
  /**
  ## A simple array */
  
  {
    double array[2] = {1.,2.};
    foreach() {
      s[] = array[0];
      s1[] = array[1];
    }
    foreach (serial)
      fprintf (stderr, "6) %g %g\n", s[], s1[]);
  }

  /**
  ## A "variable-size" array */

  {
    int nl = 2;
    double array[nl];
    array[0] = 3.; array[1] = 4.;
    foreach() {
      s[] = array[0];
      s1[] = array[1];
    }
    foreach (serial)
      fprintf (stderr, "7) %g %g\n", s[], s1[]);
  }

  /**
  ## List of scalars */

  scalar * list = {s, s1};
  foreach() {
    scalar a = list[0], b = list[1];
    a[] = 5., b[] = 6.;
  }
  foreach (serial)
    fprintf (stderr, "8) %g %g\n", s[], s1[]);
  
  foreach()
    s1[] = 7;
  
  foreach() {
    scalar a = list[0], b = list[1];
    a[] = b[];
  }
  foreach (serial)
    fprintf (stderr, "9) %g %g\n", s[], s1[]);

  /**
  ## List of vectors */

  {
    vector v[], v1[];
    vector * list = {v, v1};
    foreach() {
      vector a = list[0], b = list[1];
      foreach_dimension()
	a.x[] = 5., b.x[] = 6.;
    }
    foreach (serial)
      fprintf (stderr, "10) %g %g %g %g\n", v.x[], v.y[], v1.x[], v1.y[]);
  }

  /**
  ## For (scalar in ...) */

  foreach()
    for (scalar s in list)
      s[] = 12;
  foreach (serial)
    fprintf (stderr, "11) %g %g\n", s[], s1[]);

  /**
  ## A list with a single element */

  scalar * single = {s};
  foreach()
    for (scalar s in single)
      s[] = 12;
  
  /**
  ## An empty list */

  scalar * empty = NULL;
  foreach() {
    s1[] = 1.;
    for (scalar s in empty)
      s[] = 12;
  }
  
  /**
  ## For (s,v in list1,list2) */

  {
    vector v[], v1[];
    vector * list = {v, v1};
    scalar * list1 = {s, s1};
    foreach() {
      scalar s;
      vector v;
      for (s, v in list1, list) {
	s[] = 1;
	v.x[] = 2, v.y[] = 3;
      }
    }
    foreach (serial)
      fprintf (stderr, "12) %g %g %g %g %g %g\n", v.x[], v.y[], v1.x[], v1.y[], s[], s1[]);
  }

  /**
  ## Lists of variable length 

  This is not trivial to handle on GPUs since they do not allow for
  dynamic memory allocation. The solution is to use a different
  (compiled) static function for each list length... */

  variable_list ({s, s1});
  foreach (serial)
    fprintf (stderr, "13) %g %g\n", s[], s1[]);

  variable_list ({s1, s});
  foreach (serial)
    fprintf (stderr, "14) %g %g\n", s[], s1[]);

  scalar s2[];
  variable_list ({s, s1, s2});
  foreach (serial)
    fprintf (stderr, "15) %g %g %g\n", s[], s1[], s2[]);

  /**
  ## Functions with "inout" parameters */
  
  foreach() {
    double b = 0;
    myfunc3 (&b);
    myfunc3 (&b);
    s[] = b;
  }
  foreach (serial)
    fprintf (stderr, "16) %g\n", s[]);

  /**
  ## Constant fields and point functions */

  foreach()
    s[] = 1.;
  foreach()
    myfunc4 (point, s1, s);
  foreach (serial)
    fprintf (stderr, "17) %g %g\n", s[], s1[]);
  const scalar c[] = 2;
  foreach()
    myfunc4 (point, s1, c);
  foreach (serial)
    fprintf (stderr, "18) %g %g\n", s1[], c[]);
  
  /**
  ## Attributes */

  foreach()
    s[] = s.block + (s.nodump ? 1 : 0);
  foreach (serial)
    fprintf (stderr, "19) %g\n", s[]);
  
  /**
  ## Function pointers */
  
  double (* myfuncp) (double x) = myfunc;
  foreach()
    s[] = myfuncp (1.);
  foreach (serial)
    fprintf (stderr, "20) %g\n", s[]);

  {
    s.func = myfunc;
    s1.func = myfunc5;
    scalar * list = {s, s1};
    foreach()
      for (scalar s in list)
	s[] = s.func (2);
    foreach (serial)
      fprintf (stderr, "21) %g %g\n", s[], s1[]);
  }

  s.gradient = minmod2;
  foreach()
    if (s.gradient != NULL && s.gradient != zero)
      s1[] = s.gradient (s[-1], s[], s[1])/Delta;

  /**
  ## Functions */

  init_grid (16);
  foreach()
    s[] = myfunc (x);
  stats sf = statsf (s);
  fprintf (stderr, "22) %g %g %g\n", sf.min, sf.sum, sf.max);

  /**
  ## foreach_point() */

  init_grid (2);
  origin (-0.5, -0.5);
  foreach()
    s[] = 0.;
  for (double xp = - 0.24; xp < 0.5; xp += 0.5)
    for (double yp = - 0.24; yp < 0.5; yp += 0.5)
      foreach_point (xp, yp)
	s[] = (x + y) - (xp + yp);
  foreach (serial)
    fprintf (stderr, "23) %g %g %g\n", x, y, s[]);
  origin (0, 0);

  /**
  ## Interpolation
  
  This also tests foreach_point() and reductions. */
  
  foreach()
    s[] = x*y;
  fprintf (stderr, "24) %g %g\n",
	   interpolate (s, 0.5, 0.5, linear = false),
	   interpolate (s, 0.5, 0.5, linear = true));

  /**
  ## foreach_vertex() coordinates */

  {
    vertex scalar a[], b[];
    foreach_vertex()
      a[] = x, b[] = y;
    foreach (serial)
      fprintf (stderr, "25) %g %g\n", a[], b[]);
  }

  /**
  ## Reduction on faces */

  {
    face vector f[];
    double max = 0., sum = 0.;
    foreach_face (reduction(max:max) reduction(+:sum)) {
      f.x[] = f.x[]; // so that the loop is done on GPUs
      if (x > max) max = x;
      sum += x;
    }
    fprintf (stderr, "26) %g %g\n", sum, max);

    max = 0., sum = 0.;
    foreach_face (x, reduction(max:max) reduction(+:sum)) {
      f.x[] = f.x[]; // so that the loop is done on GPUs
      if (x > max) max = x;
      sum += x;
    }
    fprintf (stderr, "26a) %g %g\n", sum, max);
    
    max = 0., sum = 0.;
    foreach_face (y, reduction(max:max) reduction(+:sum)) {
      f.x[] = f.x[]; // so that the loop is done on GPUs
      if (x > max) max = x;
      sum += x;
    }
    fprintf (stderr, "26b) %g %g\n", sum, max);
    
    init_grid (512);
    max = 0., sum = 0.;
    foreach_face (reduction(max:max) reduction(+:sum)) {
      f.x[] = f.x[]; // so that the loop is done on GPUs
      if (x > max) max = x;
      sum += x;
    }
    fprintf (stderr, "26c) %g %g\n", sum, max);
  }

  {
    dimensions (nx = 4);
    init_grid (8);
    face vector uf[];
    foreach_face()
      uf.x[] = 1;
    double min = 100, sum = 0;
    foreach_face(reduction(min:min) reduction(+:sum)) {
      if (uf.x[] < min)
	min = uf.x[];
      sum++;
    }
    fprintf (stderr, "26e) %g %g\n", min, sum);
    dimensions (nx = 1);
    init_grid (512);
  }
  
  /**
  ## Reduction on vertices */

  {
    vertex scalar s[];
    double max = 0., sum = 0.;
    foreach_vertex (reduction(max:max) reduction(+:sum)) {
      s[] = s[]; // so that the loop is done on GPUs
      if (x > max) max = x;
      sum += x;
    }
    fprintf (stderr, "26d) %g %g\n", sum, max);
  }

  /**
  ## Staggering */

  {
    init_grid (1);
    vertex scalar a[];
    foreach() {
      double b = a[]; b = b; // fixme: a needs to be used
      v.x[] = a.d.x, v.y[] = a.d.y;
    }
    foreach (serial)
      fprintf (stderr, "27) %g %g\n", v.x[], v.y[]);
  }

  /**
  ## Multigrid */

  {
    init_grid (16);
    reset ({s}, 31.);
    foreach_level (2)
      s[] = x*y;
  
    foreach_level (2, serial)
      fprintf (stderr, "28) %g %g %g\n", x, y, s[]);

    boundary_level ({s}, 2);
  
    scalar g[];
    foreach_level_or_leaf (2)
      g[] = s[] - s[-1];
    
    foreach_level (2, serial)
      fprintf (stderr, "29) %g %g %g %g\n", x, y, g[], s[-1]);

    foreach_coarse_level (1) {
      double sum = 0.;
      foreach_child()
	sum += s[];
      s[] = sum/(1 << dimension);
    }

    foreach_level (1, serial)
      fprintf (stderr, "30) %g %g %g\n", x, y, s[]);

    foreach_level (3)
      s[] = bilinear (point, s);
    foreach_level (3, serial)
      fprintf (stderr, "31) %g %g %g\n", x, y, s[]);

    foreach()
      s[] = x*y;
    restriction ({s});
    foreach_level (2, serial)
      fprintf (stderr, "32) %g %g %g\n", x, y, s[]);
  }

  /**
  ## Boundary conditions on faces */

  {    
    init_grid (4);
    face vector uf[];
    uf.n[left] = x;
    uf.n[right] = x;
    uf.n[top] = y;
    uf.n[bottom] = y;
    foreach_face()
      uf.x[] = (x + 1)*(y + 1);
    foreach_face (x, serial)
      fprintf (stderr, "33) %g %g %g %g %g\n", x, y, uf.x[], uf.x[0,1], uf.x[0,-1]);
    foreach_face (y, serial)
      fprintf (stderr, "34) %g %g %g %g %g\n", x, y, uf.y[], uf.y[1], uf.y[-1]);    
  }

  /**
  ## Local and global variables */

  {
    scalar v[]; // same name as the global vector field 'v' above
    double vv = 2; // same name as the global double above
    foreach() {
      myfunc6 (point, v);
      myfunc7 (point, v, vv);
    }
  }

  /**
  ## Changing boundary conditions */

  {
    init_grid (1);
    reset ({s, s1}, 1);
    s[left] = neumann (0);
    myfunc8 (s);
    foreach (serial)
      assert (s[-1] == 0);
    s[left] = 1;
    myfunc8 (s);
    foreach (serial)
      assert (s[-1] == 1);
    s[left] = 2;
    scalar s1[];
    foreach()
      s1[] = s[-1];
    foreach (serial)
      assert (s1[] == 2);
  }

  /**
  ## Initialisation of `coord` */

  {
    init_grid (1);
    vector v[];
    foreach() {
      coord v0 = {1, -1}; // needs a third coordinate
      coord v1 = {0, -1, 2}, v2 = {1}; // needs a second and third coordinate
      foreach_dimension()
	v.x[] = v0.x = v1.x = v2.x;
    }
  }
  
  /**
  ## Other tests */
  
  init_grid (argc > 1 ? atoi(argv[1]) : 64);
  
  //  periodic (right);
  //  periodic (top);
  
  size (2.*pi);

  double a = 1., b = 1.;
  foreach()
    s1[] = s[] = cos(a*x)*cos(b*y);
  
  scalar p[], tmp[];
  reset ({p}, 0.);

  timer t = timer_start();
  int iter;
  for (iter = 0; iter < 40000*64/N; iter++) {

    /**
    There are two versions: the first one uses an explicit temporary
    field to avoid concurrent read/write accesses to 'p'.
    
    The second one uses concurrent read/write accesses, which usually
    does not work well on the GPU but works OK on the CPU. This is
    because in this case the convergence rate of the relaxation
    depends on the order of traversal of the cells, which will be
    different between the GPU and the CPU. */
    
#if 1    
    foreach()
      tmp[] = (p[1] + p[-1] + p[0,1] + p[0,-1] - s[]*sq(Delta))/4.;
    swap (scalar, tmp, p);
#else
    foreach()
      p[] = (p[1] + p[-1] + p[0,1] + p[0,-1] - s[]*sq(Delta))/4.;
#endif
  }

#if _GPU
  glFinish(); // make sure rendering is done on the GPU
#endif
  
  double elapsed = timer_elapsed (t);
  fprintf (stdout, "N: %d elapsed: %g speed: %g\n",
	   N, elapsed, grid->tn*iter/elapsed);

  stats sp = statsf (p);
  fprintf (stderr, "sp: %g %.5f %g\n", sp.min, fabs(sp.sum), sp.max);
  
#if 0
  foreach (serial)
    printf ("%g %g %g %g\n", x, y, p[], s1[]);
#else
  output_ppm (p, file = "p.png", n = 512, spread = -1);
  output_ppm (s1, file = "s1.png", n = 512, spread = -1);
#endif
}
