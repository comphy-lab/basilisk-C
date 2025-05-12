/**
# Upper convected derivative

<a href="http://en.wikipedia.org/wiki/Upper-convected_time_derivative">The upper convected derivative</a> is an operator to model the eulerian time derivative of a tensor $S$. It is defined by
$$
UPC(\mathbf{S}) = \mathrm{D}_t \mathbf{S}  - \mathbf{L}.\mathbf{S} - \mathbf{S}.\mathbf{L}^T
$$  
where $\mathrm{D}_t$ is the material derivative and $\mathbf{L}$ is $\mathbf{L}=(\nabla v)^T$. 
$T$ denotes the tranverse operator and the dot indicates a tensor product.

We want to discretise implicitly the linear part of a enriched UPC given by the expression, 
$$
\theta \mathbf{S} + \lambda ( \partial_t \mathbf{S} - \mathbf{L}.\mathbf{S} - \mathbf{S}.\mathbf{L}^T ) = \mathbf{R}
$$  
where $\theta$ and $\lambda$ are parameters. The above equation once discretized writes, using Einstein notation, as
$$
\left(\theta + \frac{\lambda}{\Delta t} \right) S_{ij}^{n+1} - \lambda ( L_{ik}.S^{n+1}_{kj} 
+ S^{n+1}_{ik}L_{jk}) = R_{ij} + \lambda \frac{S_{ij}^n}{\Delta t}
$$
*/

#include "poisson.h"

struct Upc {
  vector u;
  tensor S;
  tensor R;
  double dt;
  (const) scalar lambda;
  scalar theta;
  tensor L;
};

#define tensor(x) (*((tensor *)&(x)))

static void relax_upper_convected (scalar * a, scalar * b, int l, void * data)
{
  struct Upc * p = data;
  scalar theta = p->theta;
  tensor L = p->L;
  tensor S = tensor(a[0]), R = tensor(b[0]);

  foreach_level_or_leaf (l) {
    foreach_dimension ()
      S.x.x[] = (R.x.x[]+L.x.y[]*(S.x.y[]+S.y.x[]))
      /(theta[]-2.*L.x.x[]);

    S.x.y[] = (R.x.y[]+(L.y.x[]*S.x.x[]+L.x.y[]*S.y.y[]))
      /(theta[]-L.x.x[]-L.y.y[]);

    /**
    The statement below is, in reality, unnecessary since the tensor is symmetric.
    However, while this [bug](/sandbox/bugs/symm.c) is in the list... Let's be cautious! */
 
    S.y.x[] = S.x.y[]; 
  }
}

static double residual_upper_convected (scalar * a, scalar * b, 
					scalar * resl, void * data)
{
  struct Upc * p = data;
  scalar theta = p->theta;
  tensor S = tensor(a[0]), R = tensor(b[0]), res = tensor(resl[0]);
  tensor L = p->L;
  double maxres = 0.;

  foreach (reduction(max:maxres)) {
    foreach_dimension() {
      res.x.x[] = R.x.x[]- theta[]*S.x.x[] +
	(2.*L.x.x[]*S.x.x[] + L.x.y[]*(S.y.x[] + S.x.y[]));
      if (fabs (res.x.x[]) > maxres)
	maxres = fabs (res.x.x[]);
    }
  
    res.x.y[] = R.x.y[]- theta[]*S.x.y[] +
      ((L.x.x[]+L.y.y[])*S.x.y[] + L.x.y[]*S.y.y[] + L.y.x[]*S.x.x[]);
    if (fabs (res.x.y[]) > maxres)
	maxres = fabs (res.x.y[]);
    
    res.y.x[] = res.x.y[]; //Se the comment above
  }

  boundary (resl);
  return maxres;
}

mgstats upper_convected_derivative (struct Upc p)
{
  tensor L[];
  p.L = L;
  /**  I allocate a new tensor 'rhs' where I will store the right hand 
       side of the above equation although it would be spared if R 
       were already defined out of the function (p.R were not empty). 
       Note that in this case p.R (R) could be reused since it is pointing to
       an existing tensor that could be modified without any compromise. 
       I have tried but without success ! */
  tensor rhs[];
  vector u = p.u;
  tensor R = p.R;
  tensor S = p.S;
  double dt = p.dt;
    
  /** By defect, $\lambda$ and $\theta$ are assumed to be one. */ 

  if(!p.lambda.i) {
    const scalar lambda[] = 1.;
    p.lambda = lambda;
  }
  (const) scalar lambda = p.lambda;

  scalar theta[];
  if (p.theta.i) {
    (const) scalar tet = p.theta;
    foreach()
      theta[] = tet[] + lambda[]/dt;
  }
  else 
    foreach()
      theta[] = 1. + lambda[]/dt;

  p.theta = theta;
 
  foreach() {
    if (p.R.x.x.i) {
      rhs.x.y[] = R.x.y[] + lambda[]*S.x.y[]/dt;
      foreach_dimension()
	  rhs.x.x[] = R.x.x[] + lambda[]*S.x.x[]/dt;
      }
    else { // The tensor R was not passed.
      rhs.x.y[] = lambda[]*S.x.y[]/dt;
      foreach_dimension()
	rhs.x.x[] = lambda[]*S.x.x[]/dt;
    }
    rhs.y.x[] = rhs.x.y[]; //Again, see the comment above...
    foreach_dimension() {
      L.x.x[] = lambda[]*(u.x[1,0]-u.x[-1,0])/(2.*Delta);
      L.x.y[] = lambda[]*(u.x[0,1]-u.x[0,-1])/(2.*Delta);
    }
  }
  p.R = rhs;
  boundary ((scalar *) {L, theta});
  restriction ((scalar *) {L, theta});

  return mg_solve ((scalar *) {S}, (scalar *) {rhs}, 
		   residual_upper_convected, relax_upper_convected, &p);
}