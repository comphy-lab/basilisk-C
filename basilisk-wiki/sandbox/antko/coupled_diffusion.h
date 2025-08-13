/**
# Coupled reaction--diffusion equations

We here design an helper function design to solve coupled reaction-diffusion equations of the type:
\[
\left\{
\begin{array}{rl}
\theta_a \partial_t a &= \nabla\cdot (D_{aa} \nabla a) + \nabla\cdot (D_{ab} \nabla b) + \alpha_a a + \beta_a b + r_a \\
\theta_b \partial_t b &= \nabla\cdot (D_{ba} \nabla a) + \nabla\cdot (D_{bb} \nabla b) + \alpha_b a + \beta_b b + r_b \\
\end{array}
\right.
\]
which can be encountered in double diffusion problems, binary mixture cross-diffusion
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

#include "coupled_poisson.h"

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
mgstats coupled_diffusion (scalar a, scalar b, double dt,
		   face vector Daa = {{-1}},  // default 1
       face vector Dab = {{-1}},  // default 0
       face vector Dba = {{-1}},  // default 0
       face vector Dbb = {{-1}},  // default 1
		   scalar ra = {-1},         // default 0
       scalar rb = {-1},         // default 0
       scalar alphaa = {-1},  // default 0
       scalar betaa  = {-1},  // default 0
       scalar alphab = {-1},  // default 0
       scalar betab  = {-1},  // default 0
		   scalar thetaa = {-1},  // default 1
       scalar thetab = {-1})  // default 1
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar ara = automatic (ra);
  scalar arb = automatic (rb);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1. [0,1] /dt;
  (const) scalar thetaa_idt = thetaa.i >= 0 ? thetaa : idt;
  (const) scalar thetab_idt = thetab.i >= 0 ? thetab : idt;
  
  if (thetaa.i >= 0)
    foreach()
      thetaa[] *= idt[];
  if (thetab.i >= 0)
    foreach()
      thetab[] *= idt[];

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (ra.i >= 0)
    foreach()
      ara[] = thetaa_idt[]*a[] - ara[];
  else // r was not passed by the user
    foreach()
      ara[] = thetaa_idt[]*a[];
  if (rb.i >= 0)
    foreach()
      arb[] = thetab_idt[]*b[] - arb[];
  else // r was not passed by the user
    foreach()
      arb[] = thetab_idt[]*b[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda1 = thetaa_idt;
  if (alphaa.i >= 0) {
    foreach()
      alphaa[] += thetaa_idt[];
    lambda1 = alphaa;
  }
  scalar mu1 = betaa;
  scalar lambda2 = alphab;
  scalar mu2 = thetab_idt;
  if (betab.i >= 0) {
    foreach()
      betab[] += thetab_idt[];
    mu2 = betab;
  }

  /**
  Finally we solve the system. */

  return coupled_poisson (a, b, ara, arb, Daa, Dba, Dab, Dbb, lambda1, lambda2, mu1, mu2);
}