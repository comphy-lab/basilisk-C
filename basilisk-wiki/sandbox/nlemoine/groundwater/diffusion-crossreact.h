/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# System of coupled reaction-diffusion equations (cross-reaction terms only)

## Notation

Instead of solving the single equation
$$
\theta\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$

we deal with the system of $n$ equations governing the evolution of $n$ scalar fields $\left\{f_i\right\}_{1\leq i \leq n}$ :

$$
\theta_i\ \partial_t f_i \quad=\quad \nabla\cdot(D_i\nabla f_i) \quad+\quad \sum_{j=1}^n\beta_{ij} f_j \quad+\quad r_i\quad, \qquad i=1...n 
$$

Just as in the case of [simple diffusion](http://basilisk.fr/src/diffusion.h) we use a time-implicit backward Euler discretisation:

$$
\theta_i\frac{f_i^{n+1} - f_i^{n}}{dt} \quad=\quad \nabla\cdot(D_i\nabla f_i^{n+1}) \quad+\quad \sum_{j=1}^n\beta_{ij} f_j^{n+1} \quad+\quad r_i\qquad\qquad i=1...n
$$
which yields
$$
\nabla\cdot(D_i\nabla f_i^{n+1})\quad +\quad \sum_{j=1}^{n}\left(\beta_{ij}-\delta_{ij}\frac{\theta_i}{dt}\right)f_j^{n+1}\quad=\quad -\frac{\theta_i}{dt}f_i^{n}\quad -\quad r_i\qquad\qquad i=1...n 
$$

where $\delta_{ij}$ is the Kronecker symbol.

## Modified user interface for coupled Poisson–Helmholtz equations

We need to apply the generic multigrid solver to the system of Poisson–Helmholtz equations

$$
\textstyle\nabla\cdot\left(\alpha_i\nabla a_i\right)\ +\ \textstyle\sum_j\lambda_{ij} a_j \ =\ b_i\qquad\qquad i=1...n
$$

We first setup the data structure required to pass the extra parameters $\left\{\alpha_i\right\}$ and $\left\{\left\{\lambda_{ij}\right\}\right\}$.
*/
#include "poisson.h"

struct SysPoisson {
  scalar * A, * B;
  face vector * ALPHA;
  scalar ** LAMBDA ;
  double tolerance;
  int nrelax, minlevel;
  scalar * RES;
};

/**
Now we need to modify the relaxation function. We first recover the extra parameters from the data pointer.
*/
static void sysrelax (scalar * al, scalar * bl, int l, void * data)
{
  struct SysPoisson * p = (struct SysPoisson *) data;
  face vector * ALPHA = p->ALPHA;
  scalar ** LAMBDA = p->LAMBDA;
  
  /**
  We use either Jacobi (under)relaxation or we directly reuse values as soon as they are updated. For Jacobi, we need to allocate space for the new field c. Jacobi is useful mostly as it gives results which are independent of the order in which the cells are traversed, and also of the order in which evolving fields are updated in each cell in the case of coupled equations. This is not the case for the simple traversal, which means for example that results will depend on whether a tree or a multigrid is used (because cells will be traversed in a different order). The same comment applies to OpenMP or MPI parallelism. In practice however Jacobi convergence tends to be slower than simple reuse.
  */
    
  #if JACOBI
  scalar * cl = list_clone(al);
  #else
  scalar * cl = al;
  #endif  
    
  /**
  We use the face values of $\alpha$ to weight the gradients of the 5-points Laplacian operator. We get the relaxation function.
  */
  
  /* Dummy variables for iterating in lists */
  scalar a, b,c;
  scalar aj,lambdaij;
  face vector alpha;
    
  foreach_level_or_leaf (l) {
    /**
    We loop on the rows of the diffusion system
    */
    int i=0;
    for (a,b,c,alpha in al,bl,cl,ALPHA)
    {
      double n = - sq(Delta)*b[], d = 0.;
      foreach_dimension() {
        n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
        d += alpha.x[1] + alpha.x[];
      }
      int j=0;
      for(aj,lambdaij in al,LAMBDA[i])
      {
        if(j==i)
          /** $\lambda_{ii}\Delta^2$ contributes to the denominator with a negative sign*/
          d += - lambdaij[]*sq(Delta);
        else
          /** $\lambda_{ij}a_j\Delta^2$ contributes to the numerator with positive sign*/
          n += lambdaij[]*sq(Delta)*aj[];
        j++;
      } 
#if EMBED
      if (p->embed_flux) {
        double c, e = p->embed_flux (point, a, alpha, &c);
        n -= c*sq(Delta);
        d += e*sq(Delta);
      }
      if (!d)
        c[] = 0., b[] = 0.;
      else
#endif // EMBED
        c[] = n/d;
      i++;
    }
  }

  /**
  For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l) {
    for(a,c in al,cl)   
      a[] = (a[] + 2.*c[])/3.;
  }
#endif
  
#if TRASH
  for(a in al){
    scalar a1[];
    foreach_level_or_leaf (l)
      a1[] = a[];
    trash ({a});
    foreach_level_or_leaf (l)
      a[] = a1[];
  }
#endif    
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double sysresidual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  struct SysPoisson * p = (struct SysPoisson *) data;
  face vector * ALPHA = p->ALPHA;
  scalar ** LAMBDA = p->LAMBDA; 
    
  double maxres = 0.;
    
  /* Declare dummy variables for iterating in lists */
    
  scalar a,b, res;
  face vector g;
  face vector alpha;
  scalar aj,lambdaij;  
    
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
    
  face vector * gl = {NULL};
    
  for(a in al)
  {
     face vector gi = new face vector;
     foreach_dimension()
       gl = (face vector *) list_append((scalar *) gl, gi.x);
  } 

  /**
  For each row $i$ indexing an evolving variable $a_i$, compute $g_i = \alpha_i\nabla a_i$ */  
  
  foreach_face(){
    for(a,g,alpha in al,gl,ALPHA)
      g.x[] = alpha.x[]*face_gradient_x (a, 0);
  }

   /**
   Now we compute the residual $\varepsilon_i = b_i - \left(\nabla\cdot g_i\right) - \sum_j \lambda_{ij}a_j$ */
       
  foreach (reduction(max:maxres), nowarning) {
    int i=0;
    for(a,b,res in al,bl,resl)
    {
      res[] = b[];
      for(aj,lambdaij in al,LAMBDA[i])
        res[] -= lambdaij[]*aj[];
      foreach_dimension()
        res[] -= (g.x[1] - g.x[])/Delta;
#if EMBED
      if (p->embed_flux) {
        double c, e = p->embed_flux (point, a, alpha, &c);
        res[] += c - e*a[];
      }
#endif // EMBED    
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
      i++;
    }
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    int i=0;
    for(a,b,res,alpha in al,bl,resl,ALPHA)
    {
    /**
   We directly compute $\varepsilon_i = b_i - \nabla\cdot\left(\alpha_i\nabla a_i\right) - \sum_j \lambda_{ij}a_j$ */  
      
      res[] = b[];
      for(aj,lambdaij in al,LAMBDA[i])
        res[] -= lambdaij[]*aj[];
      foreach_dimension()
        res[] += (alpha.x[0]*face_gradient_x (a, 0) -
	  	  alpha.x[1]*face_gradient_x (a, 1))/Delta;  
#if EMBED
      if (p->embed_flux) {
        double c, e = p->embed_flux (point, a, alpha, &c);
        res[] += c - e*a[];
      }
#endif // EMBED
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
      i++;
    }
  }
#endif // !TREE

  fprintf(stderr,"maxres = %g\n",maxres);
  fflush(stderr);  
    
  return maxres;
}

mgstats syspoisson (scalar * al, scalar * bl, face vector * ALPHA ,
                    scalar ** LAMBDA,double tolerance = 0., int nrelax = 4,
                    int minlevel = 0, scalar * res = NULL, 
                    double (* flux) (Point, scalar, vector, double *) = NULL)
{
  /**
  We need the $\left\{\alpha_i\right\}$ and the $\left\{\left\{\lambda_{ij}\right\}\right\}$ on all levels of the grid. */
    
  restriction ((scalar *) ALPHA);
  int i=0;
  for(scalar a in al){
    restriction (LAMBDA[i]);
    i++;
  }

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct SysPoisson p = {al, bl, ALPHA, LAMBDA, tolerance, nrelax, minlevel, res };
    
  mgstats s = mg_solve (al, bl, sysresidual, sysrelax, &p,
              nrelax, res, max(1, minlevel));
    
  /**
  We restore the default. */

  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}

/**
## Modified user interface for coupled diffusion
*/

mgstats sysdiffusion (scalar * fl, double dt, face vector * Dl, scalar * rl,
                      scalar ** MAT_BETA, scalar * thetal)
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We use `rl` to store the set of r.h.s. of the Poisson--Helmholtz solver. */
 
  scalar b,f,theta,lambda;
    
  foreach() {
    for(b,f,theta in rl,fl,thetal) {
      b[] *= -1.;
      b[] -= theta[]*f[]/dt;
    }
  }
  
  /**
  FIXME: For the sake of simplicity, we use `MAT_BETA` to store the $\{\lambda_{ij}\}$'s. No verification is made as whether 
  there is any `const scalar` or `(const) scalar` in it. 
  
  We only have to modify the self-reaction term with $\lambda_{ii} = \beta_{ii}-\frac{\theta_i}{dt}$ */
    
  foreach() {
    int i=0;
    for(theta in thetal){
      lambda = MAT_BETA[i][i];
      lambda[] -= theta[]/dt;
      i++;
    }
  }
    
  /**
  Finally we solve the system. */ 
    
  return syspoisson (fl, rl, Dl, MAT_BETA);
}