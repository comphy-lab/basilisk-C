/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# System of coupled reaction-diffusion equations

## Notation

Instead of solving the single equation
$$
\theta\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$

we deal with the system of $n$ equations governing the evolution of $n$ scalar fields $\left\{f_i\right\}_{1\leq i \leq n}$ :

$$
\sum_{j=1}^n\theta_{ij}\ \partial_t f_j \quad=\quad \sum_{j=1}^n\nabla\cdot(D_{ij}\nabla f_j) \quad+\quad \sum_{j=1}^n\beta_{ij} f_j \quad+\quad r_i\quad, \qquad i=1...n 
$$

Just as in the case of [simple diffusion](http://basilisk.fr/src/diffusion.h) we use a time-implicit backward Euler discretisation:

$$
\sum_{j=1}^{n}\theta_{ij}\frac{f_j^{n+1} - f_j^{n}}{dt} \quad=\quad \sum_{j=1}^n\nabla\cdot(D_{ij}\nabla f_j^{n+1}) \quad+\quad \sum_{j=1}^n\beta_{ij} f_j^{n+1} \quad+\quad r_i\qquad\qquad i=1...n
$$
which yields
$$
\sum_{j=1}^{n}\nabla\cdot(D_{ij}\nabla f_j^{n+1})\quad +\quad \sum_{j=1}^{n}\left(\beta_{ij}-\frac{\theta_{ij}}{dt}\right)f_j^{n+1}\quad=\quad -\sum_{j=1}^{n}\frac{\theta_{ij}}{dt}f_j^{n}\quad -\quad r_i\qquad\qquad i=1...n 
$$

## Modified user interface for coupled Poisson–Helmholtz equations

We need to apply the generic multigrid solver to the system of Poisson–Helmholtz equations

$$
\textstyle\nabla\cdot\left[\sum_j\alpha_{ij}\nabla a_j\right]\ +\ \textstyle\sum_j\lambda_{ij} a_j \ =\ b_i\qquad\qquad i=1...n
$$

We first setup the data structure required to pass the extra parameters $\left\{\alpha_{ij}\right\}$ and $\left\{\lambda_{ij}\right\}$.
*/
#include "poisson.h"

struct SysPoisson {
  scalar * A, * B;
  face vector ** ALPHA;
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
  face vector ** ALPHA = p->ALPHA;
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
  In each cell, the value of the 5-points Laplacian operator can be decomposed into a part depending on the cell-centered value `a[]`, and another part depending on the neighbouring cells. For example in 2D we have at cell $(x_k,y_l)$:
  
  $$\Delta x^2\ \nabla\cdot(\alpha_{ij}\nabla a_j) \approx N_{ij}^{(k,l)} - C_{ij}^{(k,l)}a_j^{(k,l)}$$
  
  with
  
  $$N_{ij}^{(k,l)} = \alpha_{ij}^{(k+\frac{1}{2},l)}a_j^{(k+1,l)}
  + \alpha_{ij}^{(k-\frac{1}{2},l)}a_j^{(k-1,l)}
  + \alpha_{ij}^{(k,l+\frac{1}{2})}a_j^{(k,l+1)}
  + \alpha_{ij}^{(k,l-\frac{1}{2})}a_j^{(k,l-1)}$$
  
  $$C_{ij}^{(k,l)} = \alpha_{ij}^{(k+\frac{1}{2},l)}
  + \alpha_{ij}^{(k-\frac{1}{2},l)}
  + \alpha_{ij}^{(k,l+\frac{1}{2})}
  + \alpha_{ij}^{(k,l-\frac{1}{2})}$$
  
so that for each row $i$ indexing an evolving variable $a_i$ of the reaction-diffusion system, the updating equation for $a_i^{(k,l)}$ at current cell $(k,l)$ reads:  
  
  $$\left[C_{ii}^{(k,l)}-\Delta x^2\,\lambda_{ii}^{(k,l)}\right]a_i^{(k,l)} 
  = \textstyle-b_i^{(k,l)}\Delta x^2\ +\ \sum_j N_{ij}^{(k,l)}
  \ +\ \sum_{j\ne i}\Delta x^2\lambda_{ij}^{(k,l)}a_j^{(k,l)}
  \ -\ \sum_{j\ne i}C_{ij}^{(k,l)}a_j^{(k,l)}$$
  */
  
  /* Dummy variables for iterating in lists */
  scalar a, b,c;
  scalar lambda;
  face vector alpha;
    
  foreach_level_or_leaf (l) {
    /**
    We loop on the rows of the diffusion system
    */
    int row=0;
    for (b,c in bl,cl)
    {
      /**
      Now we loop on the columns of the cross-diffusion and cross-reaction coefficients matrices `ALPHA` and `LAMBDA`
      */
      double n = - sq(Delta)*b[], d = 0.;
      int col=0;
      for(a,alpha,lambda in al,ALPHA[row],LAMBDA[row]) {
        if(col==row)
            d += - lambda[]*sq(Delta);
        else
            n += lambda[]*a[]*sq(Delta);      
        foreach_dimension() {
          n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
          if(col==row) 
            d += alpha.x[1] + alpha.x[];
          else
            n += -(alpha.x[1] + alpha.x[])*a[];      
        }
        col++;
      }
      c[] = n/d;
      row++;
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
  face vector ** ALPHA = p->ALPHA;
  scalar ** LAMBDA = p->LAMBDA;
    
  double maxres = 0.;

  /* Dummy variables for iterating in lists */
  scalar a,b, res;
  face vector g;
  face vector alpha;
  scalar lambda;

#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
    
  face vector * gl = {NULL};
    
  for(a in al)
  {
     face vector gi = new face vector;
     foreach_dimension()
       gl = (face vector *) list_append((scalar *) gl, gi.x);
  } 
    
  foreach_face(){
    int row=0;
    for(g in gl){
       /**
       For each row $i$ indexing an evolving variable $a_i$, compute $g_i = \sum_j\alpha_{ij}\nabla a_j$ */
       g.x[] = 0.;
       for(a,alpha in al,ALPHA[row])
        g.x[] += alpha.x[]*face_gradient_x (a, 0);
    }
    row++;
  }
  foreach (reduction(max:maxres), nowarning) {
    int row=0;
    for(b,g,res in bl,gl,resl){
      res[] = b[];
      for(a,lambda in al,LAMBDA[row])
        res[] -= lambda[]*a[];
      foreach_dimension()
        res[] -= (g.x[1] - g.x[])/Delta;
    }
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    int row=0;
    for(b,g,res in bl,gl,resl){
      res[] = b[];
      for(a,lambda,alpha in al,LAMBDA[row],ALPHA[row]){      
        res[] -= lambda[]*a[];
        foreach_dimension()
          res[] += (alpha.x[0]*face_gradient_x (a, 0) -
                    alpha.x[1]*face_gradient_x (a, 1))/Delta;
      }
      if (fabs (res[]) > maxres)
       maxres = fabs (res[]);
    }
  }
#endif // !TREE
  return maxres; 
}

mgstats syspoisson (scalar * al, scalar * bl, face vector ** ALPHA , scalar ** LAMBDA,
                    double tolerance = 0., int nrelax = 4, int minlevel = 0, scalar * res = NULL, 
                    double (* flux) (Point, scalar, vector, double *) = NULL)
{
  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  scalar lambda;
  face vector alpha;
    
  int row=0;
  for(scalar a in al){
    for(alpha,lambda in ALPHA[row],LAMBDA[row])
      restriction ({alpha,lambda});
    row++;
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

mgstats sysdiffusion (scalar * fl, double dt, face vector ** MAT_D, scalar * rl,
                      scalar ** MAT_BETA, scalar ** MAT_THETA)
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We use `rl` to store the set of r.h.s. of the Poisson--Helmholtz solver. */
 
  scalar b,f,beta,theta,lambda;
    
  foreach() {
    int row=0;   
    for(b in rl) {
      b[] *= -1.;
      for(f,theta in fl,MAT_THETA[row])
        b[] -= theta[]*f[]/dt;
      row++;
    }
  }
  
  /**
  FIXME: For the sake of simplicity, we use `MAT_BETA` to store the $\{\lambda_{ij}\}$'s. It means that there should not be any `const scalar` in this list of lists. */
    
  int nevolv = list_len(fl);
  
  foreach() {
    for(int row=0;row<nevolv;row++){
      for(beta,theta in MAT_BETA[row],MAT_THETA[row])
        beta[] -= theta[]/dt;
    } 
  }
    
  /**
  Finally we solve the system. */

  return syspoisson (fl, rl, MAT_D, MAT_BETA);
}