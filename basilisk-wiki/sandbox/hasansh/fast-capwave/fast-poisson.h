/**
# Constant density projection method

For the constant density method we need to solve the constant
coefficient Poisson equation:
$$
\nabla ^2 p^{n+1}= \nabla \cdot \left[\left(1-
\frac{\rho_0}{\rho^{n+1}}\right)\nabla \hat{p}\right]+
\frac{\rho_0}{\Delta t} \nabla \cdot u^{*}
$$ 
After solving the Poisson equation, the velocity field will be updated
based on the new pressure field
$$
u^{n+1} = u^{*}-\Delta t \left[1/\rho_0 \nabla p^{n+1}+
\left(1/\rho^{n+1}-1/\rho_0\right)\nabla \hat{p}\right]
$$
Notations are based on the paper of [Dodd and Ferrante,
2014](http://dx.doi.org/10.1016/j.jcp.2014.05.024). */

mgstats fast_project (int method, face vector u, scalar p, scalar pn,
		      (const) face vector alpha, double rho0, double dt)
{
 
  /**
  We initialize $\hat{p}=pn$ for fastpn method and $\hat{p}=2
  p^n - p^{n-1}$ for fastpstar method. */
  
  scalar phat[];
  switch(method){
  case 1: {
    foreach()
      phat[] = p[];
    boundary ({phat});
    break;
  }
  case 2: {
    foreach() {
#if 1
      phat[] = 1.5*p[] - 0.5*pn[];
#else // this was suggested by Michael Dodd but does not seem to work
      phat[] = 2.*p[] - pn[];
#endif
      pn[] = p[];
    }
    boundary ({phat, pn});
    break;
  }
  }
  
  /**
  We calculate the term $\frac{\rho_0}{\Delta t} \nabla \cdot u^{*}$
  and store it in *div_ustar* variable. */
 
  scalar div_ustar[];
  foreach() {
    div_ustar[] = 0.;
    foreach_dimension()
      div_ustar[] += u.x[1] - u.x[];
    div_ustar[] /= dt*cm[]*Delta; 
    div_ustar[] *= rho0;
  }
  boundary ({div_ustar});

  /**
  We calculate the term $\nabla \cdot
  \left[\left(1-\frac{\rho_0}{\rho^{n+1}}\right)\nabla \hat{p}\right]$
  and store it in *div_phat* variable. */
  
  face vector f[], coef[];
  scalar div_phat[];
  foreach_face()
    coef.x[]= (1. - rho0*alpha.x[]);
#if TREE
  foreach_face()
    f.x[] = coef.x[]*(phat[] - phat[-1])/(Delta);
  boundary_flux ({f});
  foreach() {
    div_phat[] = 0.;
    foreach_dimension()
      div_phat[] += f.x[1] - f.x[];
    div_phat[] /= cm[]*Delta;
  }
#else // Cartesian
  foreach() {
    div_phat[] = 0.;
    foreach_dimension()
      div_phat[] += (coef.x[1]*(phat[1] - phat[]) -
		     coef.x[]*(phat[] - phat[-1]));
    div_phat[] /= cm[]*sq(Delta);
  }
#endif
  boundary ({div_phat});

  /**
  We construct the right hand side of Poisson equation and store it in
  variable *b*. */
  
  scalar b[];
  foreach()
    b[] = div_phat[] + div_ustar[];
  boundary ({b});
  
  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */
  
  mgstats mgp = poisson (p, b, tolerance = TOLERANCE/sq(dt));
  
  /**
  And compute $\mathbf{u}^{n+1}$ using $\mathbf{u}^*$, $p$, $p^{hat}, $ */

  foreach_face()
    u.x[] -= dt*(1./rho0*(p[] - p[-1]) +
		 (alpha.x[] - 1./rho0)*(phat[] - phat[-1]))/Delta;
  boundary ((scalar *){u});   
 
  return mgp;
}
