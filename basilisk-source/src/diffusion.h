/**
# Time-implicit discretisation of reaction--diffusion equations

We want to discretise implicitly the reaction--diffusion equation
$$
\theta\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$ 
where $\beta f + r$ is a reactive term,  $D$ is the diffusion
coefficient and $\theta$ can be a density term.

Using a time-implicit backward Euler discretisation, this can be
written
$$
\theta\frac{f^{n+1} - f^{n}}{dt} = \nabla\cdot(D\nabla f^{n+1}) + \beta
f^{n+1} + r
$$
Rearranging the terms we get
$$
\nabla\cdot(D\nabla f^{n+1}) + (\beta - \frac{\theta}{dt}) f^{n+1} =
- \frac{\theta}{dt}f^{n} - r
$$
This is a Poisson--Helmholtz problem which can be solved with a
multigrid solver. */

#include "poisson.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

trace
mgstats diffusion (scalar f, double dt,
		   face vector D = {{-1}},  // default 1
		   scalar r = {-1},         // default 0
		   scalar beta = {-1},      // default 0
		   scalar theta = {-1})     // default 1
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar ar = automatic (r);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./dt;
  (const) scalar theta_idt = theta.i >= 0 ? theta : idt;
  
  if (theta.i >= 0)
    foreach()
      theta[] *= idt[];

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (r.i >= 0)
    foreach()
      ar[] = theta_idt[]*f[] - ar[];
  else // r was not passed by the user
    foreach()
      ar[] = theta_idt[]*f[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda = theta_idt;
  if (beta.i >= 0) {
    foreach()
      beta[] += theta_idt[];
    lambda = beta;
  }

  /**
  Finally we solve the system. */

  return poisson (f, ar, D, lambda);
}
