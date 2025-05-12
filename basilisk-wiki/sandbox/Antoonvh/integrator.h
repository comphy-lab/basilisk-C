/**
# Integrate a field
*/

#include "poisson.h"

double THETA_ANGLE = 0.0;
static double residual_int (scalar * al, scalar * bl, scalar * resl, void * data) {
  double maxres = 0.;
  double st = sin(THETA_ANGLE);
  double ct = cos(THETA_ANGLE);
  int ix = ct > 0 ? 0 : -1;
  int iy = st > 0 ? 0 : -1;
  scalar a = al[0], b = bl[0], res = resl[0];
  foreach(reduction(max:maxres)) {
    double ax =  (a[1 + ix] - a[ix])/Delta; 
    double ay =  (a[0, 1 + iy] - a[0, iy])/Delta;
    res[] =  -b[] - ct*ax - st*ay ;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}
/**
with a relevant relaxation-function formulation. 
 */
static void relax_int (scalar * al, scalar * bl, int l, void * data) {
  scalar a = al[0];
  scalar b = bl[0];
  double st = sin(THETA_ANGLE);
  double ct = cos(THETA_ANGLE);
  int ix = ct > 0 ? 1 : -1;
  int iy = st > 0 ? 1 : -1;
  st = fabs(st);
  ct = fabs(ct);
  /**
     The scheme appears to benefit from a few relaxations
iterations.
   */
  for (int j = 0; j < 2; j++){
    foreach_level_or_leaf(l)
      a[] = (-b[]*Delta + ct*a[ix] + st*a[0, iy])/(ct + st);
    //a[] = (a[] + 2*(-b[]*Delta + ct*a[ix] + st*a[0, iy])/(ct + st))/3;
  }
}
/**
## User interface

Let us *hope* this procedure converges. 
 */
struct Integrate_dx {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};
  
mgstats integrate_dx (struct Integrate_dx p) {
  scalar a = p.a, b = p.b;
  return mg_solve({a}, {b}, residual_int, relax_int,
	   &p, p.nrelax, p.res, minlevel = 1);
}
