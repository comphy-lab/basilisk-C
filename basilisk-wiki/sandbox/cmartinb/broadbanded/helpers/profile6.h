/**
# Flexible solution profiling functions

It may be used to generate profiles along arbitraty interfaces. 

e.g. 

~~~literatec 
profile (all, sqrt(sq(x - X0) + sq(y - Y0)), "r_prof"); 
~~~ 

writes azimutially averaged profiles of the solution along the radial
direction (origin centered) to a file named ``r_prof`'. 

Another e.g.

~~~literatec 
profile_equi (all, x, 10, "x_prof"); 
~~~ 
 
writes 10 $y-z$ averaged solution values along the x direction
a to a file named ``x_prof`'.

It works on trees, in serial, with MPI or with openmp if compiled
*correctly*...

The functionalities rely on some helpful functions via `fractions.h`.
 */
#include "fractions.h"

/**
The first function computes the averages of the *scalar* fields in a
 `list`, along the interface that is defined by the corresponding
 volume fraction field. The function returns the average grid size
 of the cells.
 */
struct ob_av {
  scalar * list;  // List of scalar fields (mandatory)
  double * v;     // array for the values (mandatory)
  scalar c;       // Volume fractions for the interface (mandatory)
  face vector cs; // Corresponding face fraction field (optional)
};

double interface_average (struct ob_av oa){
  double area = 0, gs = 0;
  int j = 0, g = 0, m = list_len (oa.list);
  scalar c = oa.c;
  double * data = oa.v;
  for (scalar s in oa.list)
    oa.v[g++] = 0.;
  foreach(reduction(+:area) reduction(+:j)
	  reduction(+:gs) reduction (+:data[:m])) {
    if (c[] > 1e-5 && c[] < 1. - 1e-5) {
      gs += Delta;
      j++;
      coord n, p[1];
      if (oa.cs.x.i)
	n = facet_normal (point, c, oa.cs);
      else
	n = mycs (point, c);
      double alpha = plane_alpha (c[], n);
      double ar = plane_area_center(n, alpha, p)/Delta;
      foreach_dimension()
	ar *= Delta;
      area += ar;
      int g = 0;
#if dimension > 2
      coord pc = {x + Delta*p->x, y + Delta*p->y, z + Delta*p->z};
#else
      coord pc = {x + Delta*p->x, y + Delta*p->y, 0.0};
#endif
      //Point point1 = locate (pc.x, pc.y, pc.z);
      /*
      if (point1.level > 0) {
        for (scalar s in oa.list) {
          data[g++] += ar*my_interpolation(point, s, pc.x, pc.y, pc.z);
        }
      }
      */
      for (scalar s in oa.list) {
        data[g++] += ar*interpolate_linear(point, s, pc.x, pc.y, pc.z);
      }
    }
  }
  g = 0;
  if (area > 0){
    for (scalar s in oa.list)
      oa.v[g++] /= area;
  }
  if (j > 0)
    return gs/(double)j;
  else
    return L0/100.; 
}
/**
The function `profiles` uses the `interface_average()` function to
loop over the domain. It therefore requires a `phi` vertex field.
 */

struct prof {
  scalar * list;     // list of scalar field. The default is `all`
  vertex scalar phi; // The phi coordinates field (default is: phi = y)
  char * fname;      // Optional file name 
  double max;        // lower phi coordinate. The default is close to its minimum
  double min;        // upper phi coordinate. The default is close to its maximum
  double rf;         // reduction factor of query heights. The default is 1
  FILE * fp;         // File pointer, if `fname` is not provided. The default is `stdout`
  int n;             // Number of phi-values to query for equidistant profiles
  char * func_str;   // An identifier string (default is `phi`)
};

void profiles (struct prof p) {
  /** 
We start with setting a bunch of default settings. Such as the list
includes all scalars, except `phi`. 
  */
  if(!p.list)
    p.list = all;
  scalar * list = NULL ;
  for (scalar s in p.list){
    if (strcmp ("phi", s.name))
	list = list_append (list, s);
  }
  vertex scalar phi;
  if (p.phi.i)
    phi = p.phi;
  /**
For historic reasons, we make a vertical profile is `phi` is not
provided.
   */
  else {
    phi = new_vertex_scalar ("phi");
    foreach_vertex ()
      phi[] = y;
    p.func_str = "y";
  }
  /**
The `min` and `max` coordinates are the extreme values, taken from
the `phi` field with the addional `maxdepth` trick of Vincent Heusinkveld.
   */
  if (!p.min){
    double mn = HUGE;
    foreach_vertex(reduction(min:mn)){
      if (phi[] < mn)
	mn = phi[];
    }
    p.min = mn + L0/(double)(1 << (grid->maxdepth + 1));
  }
  if (!p.max){
    double mx = -HUGE;
    foreach_vertex(reduction(max:mx)){
      if (phi[] > mx)
	mx = phi[];
    }
    p.max = mx - L0/(double)(1 << (grid->maxdepth + 1));
  }
  /**
     And by default there is no reduction factor and we print to the
     `stdout` file pointer.
   */
  if (!p.rf)
    p.rf = 1;
  if (!p.fname && !p.fp)
    p.fp = stdout;
  double dzn;
  if (p.n)
    dzn = (p.max - p.min) / ((double)p.n - 0.99999999);
  if (!p.func_str)
    p.func_str = "phi";
  FILE * fp = p.fp;
  char * file = p.fname;
  if (pid()==0){ 
    if (file && (fp = fopen (file, "w")) == NULL) {
      perror (file);
      exit (1);
    }
    assert (fp);
    fprintf(fp,"#%s\t", p.func_str);
    for(scalar s in list)
      fprintf(fp,"%s\t",s.name);
    fprintf(fp,"\n");
  }
  int len = list_len(list);
  boundary(list);
  double phi_p = p.min;
  scalar s[];
  face vector sf[];
  while (phi_p < p.max){
    double aver[len];
    fractions (phi, s, sf, val = phi_p);
    double dphi = interface_average (list, aver, s, sf);
    if (pid() == 0){
      int k = 0;
      for (scalar s in list){
	if (k == 0)
	  //fprintf (fp, "%g\t%g", phi_p, aver[k]);
	  fprintf (fp, "%8E\t%8E", phi_p, aver[k]);
	if (k > 0)
	  //fprintf(fp,"\t%g",aver[k]);
	  fprintf(fp,"\t%8E", aver[k]);
	k++;
	if (k == len)
	  fprintf(fp,"\n");
      }
    }
    if (p.n)
      dphi = dzn/p.rf;
    /**
       Note that the units of `phi` and `dphi` may be incompatible.
     */
    phi_p += p.rf*dphi;
  }
  if (pid()==0){
    fflush(fp);
    if (fp != stdout)
      fclose(fp);
  }
  free (list);
}
/**
## Two addional user-interface functions

Because the `phi` field may be inconvienient, we make some
user-friendly defenitions. Note that I do not know how to do optional
arguments, so the user *must* provide all arguments. It relies on the
default settings of the `profiles()` function above.

For grid-adaptive profiles use `profile()`:
 */

#define profile(list, fun, fname) do {		\
    vertex scalar phi[];			\
    foreach_vertex()					\
      phi[] = fun;					\
    char * str = #fun;					\
    profiles (list, phi, fname, func_str = str);	\
  } while(0);

/**
And for equidistant profiles use `profile_equi()`:
*/
#define profile_equi(list, fun, nr, fname) do {	\
    vertex scalar phi[];			\
    foreach_vertex()				\
      phi[] = fun;					\
    char * str = #fun;					\
    profiles (list, phi, fname, n = nr, func_str = str);	\
  } while(0);
