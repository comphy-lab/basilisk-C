/**
# Combined Implicit-Expicit Runge-Kutta Scheme

One can use Implicit-Explicit (IMEX) splitting to solve DE's.

*/
#include "run.h"
#include "utils.h"

/**
The tendency of the explicit terms need to be computed in this
function:
 */

extern void explicit_tendency (scalar s, scalar ds);

/**
The user must define a method to solve the implicit problem for `ds`$_{imp}
=g(s)$:

$$s_{new} = s_{old} + dt A_{i,i} g(s_{new}) $$

The function must *both* update $s \leftarrow s_{new}$ and compute the
tendency (see example/tests).
*/

extern void implicit_tendency (scalar s, scalar ds,
			       double dtAii);

/**
## Some Butcher arrays

The solver can solve for every Explicit + Diagonally implicit RK
method. We provide two:
 */

#if IMEX_222 
#define RKS 2
// Explicit
double A[RKS][RKS] = {{0, 0},
		      {1./2., 0}};
double b[RKS] = {0, 1};
// Implicit
double Ai[RKS][RKS] = {{0, 0},
		       {0, 1./2.}};
double bi[RKS] = {0, 1};
#endif 

#if IMEX_333
#define RKS 3
#define SQRT3  1.73205080756887729353
#define GAMMA ((3 + SQRT3)/6.)
// Explicit
double A[RKS][RKS] = {{0, 0, 0},
		      {GAMMA, 0, 0},
		      {GAMMA - 1., 2.*(1. - GAMMA), 0}};
double b[RKS] = {0, 1./2., 1./2.};
// Explicit
double Ai[RKS][RKS] = {{0, 0, 0},
		       {0, GAMMA, 0},
		       {0, 1 - 2*GAMMA, GAMMA}};
double bi[RKS] = {0, 1./2., 1./2.};
#endif
/**
Peform two tests on the arrays
 */
bool test_method() {
  for (int i = 0; i < RKS; i++) {
    for (int j = 0; j < RKS; j++) {
      if (i >= j && A[j][i]) {
	fprintf (stderr, "The A matrix is not an explicit method\n");
	return false;
      }
      if (i > j && Ai[j][i]) {
	fprintf (stderr, "The Ai matrix is not an diagonally implicit method\n");
	return false;
      }
    }
  }
  return true;
}
/**
2 $\times$ `RKS` fiels are allocated to store the tendencies.
 */
scalar * dsl = NULL, * dsiml = NULL;
void scalars_for_imex (void) {
  for (int i = 0; i < RKS; i++) {
    dsl = list_add (dsl, new_scalar ("dsi"));
    dsiml = list_add (dsiml, new_scalar ("dsim"));
  }
}

void cleanup_scalars_for_imex (void) {
  delete (dsl); free (dsl); dsl = NULL;
  delete (dsiml); free (dsiml); dsiml = NULL;
}

event defaults (i = 0) {
  if (!test_method())
    return 1;
  scalars_for_imex();
}

event clean (t = end) {
  cleanup_scalars_for_imex();
}

event time_step(i++);
/**
The IMEX-timestep event is implemented below
 */
extern scalar * sl;
event IMEX_advance (i++, last) {
  scalar s_temp[];
  scalar s = sl[0];
  for (int i = 0; i < RKS; i++) {
    // explicit terms;
    foreach()
      s_temp[] = s[];
    for (int j = 0; j < i; j++) {
      if (A[i][j]) {
	scalar dsj = dsl[j];
	double Aji = A[i][j];
	foreach()
	  s_temp[] += dt*Aji*dsj[];
      }
      if (Ai[i][j]) {
     
	scalar dsj = dsiml[j];
	double Aiji = Ai[i][j];
	foreach()
	  s_temp[] += dt*Aiji*dsj[];
      }
    }
    // Is this an implicit stage? 
    if (Ai[i][i]) {
      // Find implicit tendency
      scalar dsimi = dsiml[i];
      implicit_tendency (s_temp, dsimi, dt*Ai[i][i]);
    }
    // Compute the explicit tendency
    scalar dsi = dsl[i];
    explicit_tendency (s_temp, dsi);
  }
  foreach() 
    for (int i = 0; i < RKS; i++) {
      if (b[i]) {
	scalar dsi = dsl[i];
	s[] += dt*b[i]*dsi[];
      }
      if (bi[i]) {
	scalar dsimi = dsiml[i];
	s[] += dt*bi[i]*dsimi[];
      }
    }
}


/**
## Test
* [A test](imex.c)

## See also
* [Implicit Gauss-Legendre Runge-Kutta methods](GLrk.h)
*/
