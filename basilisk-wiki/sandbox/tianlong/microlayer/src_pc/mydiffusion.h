#ifndef _MY_DIFFUSION_H
#define _MY_DIFFUSION_H
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
#if USE_CONJUGATE_HEAT
#include "poissonconjugate.h"
#endif

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */
double TOL_DIFFUSION = 1e-3;
struct Diffusion {
  // mandatory
  scalar f;
  double dt;
  // optional
  face vector D;  // default 1
  scalar r, beta; // default 0
  scalar theta;   // default 1
  scalar *res;
};

trace
mgstats diffusion (struct Diffusion p)
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar f = p.f, r = automatic (p.r);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./p.dt;
  (const) scalar theta_idt = p.theta.i ? p.theta : idt;
  
  if (p.theta.i) {
    scalar theta_idt = p.theta;
    foreach()
      theta_idt[] *= idt[];
  }

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.r.i)
    foreach()
      r[] = theta_idt[]*f[] - r[];
  else // r was not passed by the user
    foreach()
      r[] = theta_idt[]*f[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda = theta_idt;
  if (p.beta.i) {
    scalar beta = p.beta;
    foreach()
      beta[] += theta_idt[];
    lambda = beta;
  }

  /**
  Finally we solve the system. */

  return poisson (f, r, p.D, lambda, tolerance = TOL_DIFFUSION, res = p.res);
}

#if USE_CONJUGATE_HEAT
struct myDiffusion
{
  // mandatory
  scalar f1, f2;
  double dt;
  // optional
  face vector D1, D2;          // default 1
  scalar r1, r2, beta1, beta2; // default 0
  scalar theta1, theta2;       // default 1
};

trace
mgstats mydiffusion (struct myDiffusion p)
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar fl = p.f1, rl = automatic (p.r1);
  scalar fg = p.f2, rg = automatic (p.r2);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./p.dt;
  (const) scalar thetal_idt = p.theta1.i ? p.theta1 : idt;
  (const) scalar thetag_idt = p.theta2.i ? p.theta2 : idt;
  
  if (p.theta1.i) {
    scalar thetal_idt = p.theta1;
    foreach()
      thetal_idt[] *= idt[];
  }

  if (p.theta2.i) {
    scalar thetag_idt = p.theta2;
    foreach()
      thetag_idt[] *= idt[];
  }

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.r1.i)
    foreach()
      rl[] = thetal_idt[]*fl[] - rl[];
  else // r was not passed by the user
    foreach()
      rl[] = thetal_idt[]*fl[];

  if (p.r2.i)
    foreach()
      rg[] = thetag_idt[]*fg[] - rg[];
  else // r was not passed by the user
    foreach()
      rg[] = thetag_idt[]*fg[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambdal = thetal_idt;
  scalar lambdag = thetag_idt;

  if (p.beta1.i) {
    scalar betal = p.beta1;
    foreach()
      betal[] += thetal_idt[];
    lambdal = betal;
  }

  if (p.beta2.i) {
    scalar betag = p.beta2;
    foreach()
      betag[] += thetag_idt[];
    lambdag = betag;
  }

  scalar *a = {fl, fg};
  scalar *b = {rl, rg};
  face vector *alpha = {p.D1, p.D2};
  scalar *lambda = {lambdal, lambdag};
  mgstats solution = mypoisson(a, b, alpha, lambda, tolerance = TOL_DIFFUSION);

  return solution;
}


struct myDiffusionConjugate
{
  // mandatory
  scalar fl, fg, fs;
  double dt;
  // optional
  face vector Dl, Dg, Ds;          // default 1
  scalar rl, rg, rs, betal, betag, betas; // default 0
  scalar thetal, thetag, thetas;       // default 1
  scalar resl, resg, ress;
};

trace
mgstats mydiffusionConjugate (struct myDiffusionConjugate p)
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar fl = p.fl, rl = automatic (p.rl);
  scalar fg = p.fg, rg = automatic (p.rg);
  scalar fs = p.fs, rs = automatic (p.rs);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./p.dt;
  (const) scalar thetal_idt = p.thetal.i ? p.thetal : idt;
  (const) scalar thetag_idt = p.thetag.i ? p.thetag : idt;
  (const) scalar thetas_idt = p.thetas.i ? p.thetas : idt;
  
  if (p.thetal.i) {
    scalar thetal_idt = p.thetal;
    foreach()
      thetal_idt[] *= idt[];
  }

  if (p.thetag.i) {
    scalar thetag_idt = p.thetag;
    foreach()
      thetag_idt[] *= idt[];
  }

  if (p.thetas.i) {
    scalar thetas_idt = p.thetas;
    foreach()
      thetas_idt[] *= idt[];
  }

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.rl.i)
    foreach()
      rl[] = thetal_idt[]*fl[] - rl[];
  else // r was not passed by the user
    foreach()
      rl[] = thetal_idt[]*fl[];

  if (p.rg.i)
    foreach()
      rg[] = thetag_idt[]*fg[] - rg[];
  else // r was not passed by the user
    foreach()
      rg[] = thetag_idt[]*fg[];

  if (p.rs.i)
    foreach()
      rs[] = thetas_idt[]*fs[] - rs[];
  else // r was not passed by the user
    foreach()
      rs[] = thetas_idt[]*fs[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambdal = thetal_idt;
  scalar lambdag = thetag_idt;
  scalar lambdas = thetas_idt;

  if (p.betal.i) {
    scalar betal = p.betal;
    foreach()
      betal[] += thetal_idt[];
    lambdal = betal;
  }

  if (p.betag.i) {
    scalar betag = p.betag;
    foreach()
      betag[] += thetag_idt[];
    lambdag = betag;
  }

  if (p.betas.i) {
    scalar betas = p.betas;
    foreach()
      betas[] += thetas_idt[];
    lambdas = betas;
  }

  scalar *a = {fl, fg, fs};
  scalar *b = {rl, rg, rs};
  face vector *alpha = {p.Dl, p.Dg, p.Ds};
  scalar *lambda = {lambdal, lambdag, lambdas};
  scalar *res = {p.resl, p.resg, p.ress};
  mgstats solution = mypoisson(a, b, alpha, lambda, tolerance = TOL_DIFFUSION, res = res);

  return solution;
}
#endif //USE_CONJUGATE_HEAT

#endif //MY_DIFFUSION_H