/**
# Evaluate a line integral 

We want to find $a$, such that,

$$\frac{\partial a}{\partial n} = -b$$

With suitable boundary conditions, this corresponds to, 

$$a = \int b \mathrm{d}n$$

with $n$ some coordinate line under an angle $\theta$ (`THETA_ANGLE`) with the $x$ axis.
For this purpose a multigrid solver is used,
*/

#include "grid/multigrid.h"
#include "poisson.h"
/**
A residual function for our problem, 
 */
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
struct Integrate_dn {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};
  
mgstats integrate_dn (struct Integrate_dn p) {
  scalar a = p.a, b = p.b;
  return mg_solve({a}, {b}, residual_int, relax_int,
	   &p, p.nrelax, p.res, minlevel = 1);
}

/**
## Give it a go.
 */

#include "utils.h"
scalar a[], b[];
a[right] = dirichlet (0);
a[top] = dirichlet (0);
a[left] = dirichlet (0);
a[bottom] = dirichlet (0);

int main() {
  X0 = Y0 = -L0/2.;
  init_grid (1 << 9);
  foreach()
    b[] = exp(-(sq(x*10) + sq(y*10)));
  output_ppm (b, file = "b.png", n = 256, min = 0, max = 1);
  restriction({b});
  boundary ({b});
  for (double THETA = 0.; THETA < 3*pi; THETA += 0.1){
    THETA_ANGLE = THETA;
    printf("%d\n", integrate_dn (a, b).i);
    output_ppm (a, file = "a.mp4", n = 256, min = 0, max = sqrt(pi)/10.);
  }
}

/**
# Results

The source term ($b$) looks like this:
   
![The $b$ field](my_first_MG_solver/b.png)

And the solution ($a$):

![The $a$ field](my_first_MG_solver/a.mp4)

Well done multigrid solver! Much better than [this](smoke.c) $\mathcal{O}(N^2)$ approach! 
 */
