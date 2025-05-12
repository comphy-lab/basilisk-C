/**
# Vector potential solver

The edge-vector potential ($\vec A$) for a face vector field $\vec
{R}$. If $ \vec R$ is not Solenoidal, the vector potential will
"simply" give the divergence-free part. Because edges do not really
exist in Basilisk, it only works for full tree grids.

The vector potential ($\vec{A}$) is found iteratively, by solving
$$\nabla_D \times \vec{F} = \nabla_D \times \nabla \times \vec A,$$
where "$\nabla_D \times$" is the approximate curl (from faces to
edge).

Note: The MG-cycle is not designed to work for non-cell centered
solutions. It would be nice if we could give it the relevant
cache(level)s and prolongation method.
 */

#include "poisson.h"
#include "utils.h"

double residual_v (scalar * Al, scalar * omegal,
		   scalar * resl, void * data) {
  vector A, omega, res;
  int i = 0;
  foreach_dimension() {
    A.x = Al[i];
    omega.x = omegal[i];
    res.x = resl[i++];
  }
  double max_res = 0;
  foreach(reduction(max:max_res)) {
    foreach_dimension() {
      double F1y = (A.z[] + A.x[0,0,1] - A.z[1] - A.x[])/Delta; 
      double F2y = (A.z[] - A.x[-1,0,1] - A.z[-1] + A.x[-1])/Delta; 
      double F1x = (A.z[] + A.y[0,0,1] - A.z[0,1] - A.y[])/Delta; 
      double F2x = (A.z[] - A.y[0,-1,1] - A.z[0,-1] + A.y[0,-1])/Delta;
      res.z[] = omega.z[] - (F1x + F2x + F1y + F2y)/Delta;
    }
    double resi = max(max(fabs(res.x[]), fabs(res.y[])), fabs(res.z[]));
    if (resi > max_res)
      max_res = resi;
  }
  return max_res;
}

double under_relax = 0.6;

void relax_v (scalar * Al, scalar * resl,
		int l, void * data) {
  vector A, res;
  int i = 0;
  foreach_dimension() {
    A.x = Al[i];
    res.x = resl[i++];
  }
  foreach_level_or_leaf(l) {
    foreach_dimension()
      A.z[] = under_relax*res.z[]*sq(Delta)/4.;
  }
}

void vector_potential (vector R, vector A, bool Coulomb_gauge = false) {
  vector omg[];
  foreach() {
    foreach_dimension() 
      omg.z[] = (R.y[] - R.x[] - R.y[-1] + R.x[0,-1])/Delta;
  }
  mg_solve ((scalar*){A}, (scalar*){omg}, residual_v, relax_v);

  stats f_x = statsf(A.x);
  stats f_y = statsf(A.y);
  stats f_z = statsf(A.z);
  foreach() {
    foreach_dimension() {
      A.x[] -= f_x.sum/(cube(L0));
    }
  }
  scalar div[], p[];
  if (Coulomb_gauge) { // vertex divergence
    foreach() {
      div[] = 0;
      p[] = 0;
      foreach_dimension()
	div[] += (A.x[] - A.x[-1])/Delta;
    }
    poisson(p, div);
    foreach() {
      foreach_dimension()
	A.x[] -= (p[1] - p[])/Delta;
    }
  }
}
/**
## Test

* [Helmholtz-Hodge decomposition of a noisy 3D vectorfield](HHD.c)
*/