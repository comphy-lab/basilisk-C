/**
# A solver for coupled Poisson equations

We want to the following set of coupled Poisson--Helmholtz equations:
\[
\left\{
\begin{array}{rl}
\nabla\cdot (\alpha_1 \nabla a) + \nabla\cdot (\beta_1 \nabla b) + \lambda_1 a + \mu_1 b &= c \\
\nabla\cdot (\alpha_2 \nabla a) + \nabla\cdot (\beta_2 \nabla b) + \lambda_2 a + \mu_2 b &= d
\end{array}
\right.
\]
To do so, we will build on the generic multigrid solver [poisson.h](poisson.h).
*/

#include "poisson.h"

/**
## Helper function to solve a 2x2 system
At each point we will need to solve a 2x2 system of equations. The following helper function performs this operation.
 */
static void solve_2x2 (const double matrix[2][2], const double rhs[2], double result[2])
{
  double determinant = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
  result[0] = (rhs[0]*matrix[1][1] - rhs[1]*matrix[0][1]) / determinant;
  result[1] = (matrix[0][0]*rhs[1] - matrix[1][0]*rhs[0]) / determinant;
}

struct CoupledPoisson {
  scalar a, b, c, d;
  (const) face vector alpha1, alpha2, beta1, beta2;
  (const) scalar lambda1, lambda2, mu1, mu2;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

/**
## Relaxation function
We now write the relaxation function using the generic call structure for `relax()`.
 */

static void relax_coupled_poisson (scalar * unknowns, scalar * rhs, int l, void * data)
{
  scalar a = unknowns[0], b = unknowns[1];
  scalar c = rhs[0], d = rhs[1];
  struct CoupledPoisson * cp = (struct CoupledPoisson *) data;
  (const) face vector alpha1 = cp->alpha1, alpha2 = cp->alpha2;
  (const) face vector beta1 = cp->beta1, beta2 = cp->beta2;
  (const) scalar lambda1 = cp->lambda1, lambda2 = cp->lambda2;
  (const) scalar mu1 = cp->mu1, mu2 = cp->mu2;
  double matrixD[2][2], vectorN[2], sol[2];

  /**
  In a large majority of cases, we will make use of Gauss-Seidel relaxation. However, we also provide an option for Jacobi relaxation. 
  */
  
#if JACOBI
  scalar da[], db[];
#else
  scalar da = a, db = b;
#endif

  /**
  We keep the same loop declaration introduced in [poisson.h](poisson.h) 
  allowing for red/black Gauss-Seidel relaxation. */
  
#if GAUSS_SEIDEL || _GPU
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l, nowarning)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
  foreach_level_or_leaf (l, nowarning)
#endif
  {
    matrixD[0][0] = -lambda1[]*sq(Delta);
    matrixD[0][1] = -mu1[]*sq(Delta);
    matrixD[1][0] = -lambda2[]*sq(Delta);
    matrixD[1][1] = -mu2[]*sq(Delta);
    vectorN[0] = -c[]*sq(Delta);
    vectorN[1] = -d[]*sq(Delta);
    foreach_dimension() {
      matrixD[0][0] += alpha1.x[1] + alpha1.x[];
      matrixD[0][1] += beta1.x[1] + beta1.x[];
      matrixD[1][0] += alpha2.x[1] + alpha2.x[];
      matrixD[1][1] += beta2.x[1] + beta2.x[];
      vectorN[0] += alpha1.x[1]*a[1] + alpha1.x[]*a[-1] +
                    beta1.x[1]*b[1] + beta1.x[]*b[-1];
      vectorN[1] += alpha2.x[1]*a[1] + alpha2.x[]*a[-1] +
                    beta2.x[1]*b[1] + beta2.x[]*b[-1];
    }
    solve_2x2 (matrixD, vectorN, sol);
    da[] = sol[0];
    db[] = sol[1];
  }

  /**
  For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l) {
    a[] = (a[] + 2.*da[])/3.;
    b[] = (b[] + 2.*db[])/3.;
  }
#endif
  
#if TRASH
  scalar a1[], b1[];
  foreach_level_or_leaf (l) {
    a1[] = a[];
    b1[] = b[];
  }
  trash ({a, b});
  foreach_level_or_leaf (l) {
    a[] = a1[];
    b[] = b1[];
  }
#endif
}

/**
We design the residual function by making use of the recommendation for tree
meshes noted in (poisson.h)[poisson.h] */

static double residual_coupled_poisson (scalar * variables, scalar * rhsl, scalar * resl, void * data)
{
  scalar a = variables[0], b = variables[1];
  scalar c = rhsl[0], d = rhsl[1];
  scalar res1 = resl[0], res2 = resl[1];
  struct CoupledPoisson * cp = (struct CoupledPoisson *) data;
  (const) face vector alpha1 = cp->alpha1, alpha2 = cp->alpha2;
  (const) face vector beta1 = cp->beta1, beta2 = cp->beta2;
  (const) scalar lambda1 = cp->lambda1, lambda2 = cp->lambda2;
  (const) scalar mu1 = cp->mu1, mu2 = cp->mu2;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector ga1[], gb1[], ga2[], gb2[];
  foreach_face() {
    ga1.x[] = alpha1.x[]*face_gradient_x (a, 0);
    gb1.x[] = beta1.x[] *face_gradient_x (b, 0);
    ga2.x[] = alpha2.x[]*face_gradient_x (a, 0);
    gb2.x[] = beta2.x[] *face_gradient_x (b, 0);
  }
  foreach (reduction(max:maxres), nowarning) {
    res1[] = c[] - lambda1[]*a[] - mu1[]*b[];
    res2[] = d[] - lambda2[]*a[] - mu2[]*b[];
    foreach_dimension() {
      res1[] -= (ga1.x[1] - ga1.x[])/Delta +
                (gb1.x[1] - gb1.x[])/Delta;
      res2[] -= (ga2.x[1] - ga2.x[])/Delta +
                (gb2.x[1] - gb2.x[])/Delta;
    }
    if ((fabs (res1[]) + fabs (res2[])) > maxres)
      maxres = (fabs (res1[]) + fabs (res2[]));
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    res1[] = c[] - lambda1[]*a[] - mu1[]*b[];
    res2[] = d[] - lambda2[]*a[] - mu2[]*b[];
    foreach_dimension() {
      res1[] += (alpha1.x[0]*face_gradient_x (a, 0) -
		             alpha1.x[1]*face_gradient_x (a, 1))/Delta +
                (beta1.x[0]*face_gradient_x (b, 0) -
                 beta1.x[1]*face_gradient_x (b, 1))/Delta; 
      res2[] += (alpha2.x[0]*face_gradient_x (a, 0) -
                 alpha2.x[1]*face_gradient_x (a, 1))/Delta +
                (beta2.x[0]*face_gradient_x (b, 0) -
                 beta2.x[1]*face_gradient_x (b, 1))/Delta;
    }
    if ((fabs (res1[]) + fabs (res2[])) > maxres)
      maxres = (fabs (res1[]) + fabs (res2[]));
  }
#endif // !TREE
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

trace
mgstats coupled_poisson (scalar a, scalar b, scalar c, scalar d,
		 (const) face vector alpha1 = {{-1}},
     (const) face vector alpha2 = {{-1}},
     (const) face vector beta1 = {{-1}},
     (const) face vector beta2 = {{-1}},
		 (const) scalar lambda1 = {-1},
     (const) scalar lambda2 = {-1},
     (const) scalar mu1 = {-1},
     (const) scalar mu2 = {1},
		 double tolerance = 0.,
		 int nrelax = 4,
		 int minlevel = 0,
		 scalar * res = NULL)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (alpha1.x.i < 0)
    alpha1[] = {1.,1.,1.};
  if (beta1.x.i < 0)
    alpha1[] = {0.,0.,0.};
  if (alpha2.x.i < 0)
    alpha2[] = {0.,0.,0.};
  if (beta2.x.i < 0)
    beta2[] = {1.,1.,1.};
  if (lambda1.i < 0)
    lambda1[] = 0.;
  if (lambda2.i < 0)
    lambda2[] = 0.;
  if (mu1.i < 0)
    mu1[] = 0.;
  if (mu2.i < 0)
    mu2[] = 0.;
  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  restriction ({alpha1,alpha2,beta1,beta2,lambda1,lambda2,mu1,mu2});

  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct CoupledPoisson cp = {a, b, c, d, alpha1, alpha2, beta1, beta2, lambda1, lambda2, mu1, mu2, tolerance, nrelax, minlevel, res };

  mgstats s = mg_solve ({a,b}, {c,d}, residual_coupled_poisson, relax_coupled_poisson, &cp,
			nrelax, res, max(1, minlevel));

  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}

