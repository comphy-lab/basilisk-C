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

#include "mypoisson.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

struct Diffusion {
  // mandatory
  scalar f;
  double dt;
  // optional
  face vector D;  // default 1
  scalar r, beta; // default 0
  scalar theta;   // default 1
};

trace
mgstats diffusion (scalar f, double dt,
		               face vector D = {{-1}},  // default 1
		               scalar r = {-1},         // default 0
		               scalar beta = {-1},      // default 0
		               scalar theta = {-1})     // default 1
{
  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  scalar ar = automatic (r);

  // theta_idt = -theta/dt (or -1/dt if theta not provided)
  const scalar idt[] = -1./dt;
  //(const) scalar theta_idt = theta.i >= 0 ? theta : idt;
  scalar theta_idt[];

  if (theta.i >= 0) {
    foreach()
      theta_idt[] = theta[] * idt[];
  } else {
    foreach()
      theta_idt[] = idt[];
  }
  boundary ({theta_idt});

  if (r.i >= 0)
    foreach()
      ar[] = theta_idt[]*f[] - ar[];
  else
    foreach()
      ar[] = theta_idt[]*f[];
#if EMBED
  foreach()
    ar[] *= cs[];
#endif 
  boundary ({ar});

  scalar lambda[];
  foreach()
    lambda[] = theta_idt[];
  if (beta.i >= 0) {
    foreach()
      lambda[] += beta[];
  }

#if EMBED
  foreach()
    lambda[] *= cs[];
#endif 
  boundary ({lambda});

  return poisson (f, ar, D, lambda);
}
