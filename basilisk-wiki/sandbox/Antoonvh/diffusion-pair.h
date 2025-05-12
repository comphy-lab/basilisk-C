/**
# Time-implicit discretisation of reaction--diffusion equations

This is an adaptation of [`diffusion.h`](/src/diffusion.h) for pairwise duffusion and flux partitioning

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

#include "poisson-pair.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

struct Diffusion_pair {
  // mandatory
  scalar f1;
  scalar f2;
  double dt;
  // optional
  face vector D1;  // default 1
  face vector D2;  // default 1

  scalar r1, beta1; // default 0
  scalar r2, beta2; // default 0

  scalar theta1;   // default 1
  scalar theta2;   // default 1

};

trace
mgstats diffusion_pair (struct Diffusion_pair p)
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */
  
  scalar f1 = p.f1, r1 = automatic (p.r1);
  scalar f2 = p.f2, r2 = automatic (p.r2);
  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./p.dt;
  (const) scalar theta_idt1 = p.theta1.i ? p.theta1 : idt;
  (const) scalar theta_idt2 = p.theta2.i ? p.theta2 : idt;
  if (p.theta1.i) {
    scalar theta_idt1 = p.theta1;
    foreach()
      theta_idt1[] *= idt[];
  }
  if (p.theta2.i) {
    scalar theta_idt2 = p.theta2;
    foreach()
      theta_idt2[] *= idt[];
  }

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.r1.i)
    foreach()
      r1[] = theta_idt1[]*f1[] - r1[];
  else // r was not passed by the user
    foreach()
      r1[] = theta_idt1[]*f1[];
  if (p.r2.i)
    foreach()
      r2[] = theta_idt2[]*f2[] - r2[];
  else // r was not passed by the user
    foreach()
      r2[] = theta_idt2[]*f2[];
  
  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda1 = theta_idt1;
  scalar lambda2 = theta_idt2;
  if (p.beta1.i) {
    scalar beta1 = p.beta1;
    foreach()
      beta1[] += theta_idt1[];
    lambda1 = beta1;
  }
  if (p.beta2.i) {
    scalar beta2 = p.beta2;
    foreach()
      beta2[] += theta_idt2[];
    lambda2 = beta2;
  }
  /**
  Finally we solve the system. */

  return poisson_pair ({f1, f2}, {r1, r2},
		       p.D1, p.D2, lambda1, lambda2);
}
