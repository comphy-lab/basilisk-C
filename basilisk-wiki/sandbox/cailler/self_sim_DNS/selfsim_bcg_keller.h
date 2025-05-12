/**
# SELF-SIMILAR Bell-Colella-Glaz advection scheme (Keller \& Miksis problem)

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Colella-Glaz, 1989](references.bib#bell89), **but modified to take 
into account specific source terms from non-divergence-free advection 
velocities**. 
Given a centered scalar field *f*, a face vector field *uf* (possibly weighted by a
face metric), a timestep *dt*, the number of dimension of the simulation *Nd* 
and a source term field *src*, it fills the face vector field *flux* with 
the components of the advection fluxes of *f*. */

void selfsim_tracer_fluxes (scalar f,
		    face vector uf, 
		    face vector flux,
		    double dt,
        double Nd, 
		    (const) scalar src)
{

  /**
  We first compute the cell-centered gradient of *f* in a locally-allocated
  vector field. */
  
  vector g[];
  gradients ({f}, {g});

  /**
  For each face, the flux is composed of two parts... */

  foreach_face() {

    /**
    A normal component... (Note that we cheat a bit here, `un` should
    strictly be `dt*(uf.x[i] + uf.x[i+1])/((fm.x[] +
    fm.x[i+1])*Delta)` but this causes trouble with boundary
    conditions (when using narrow '1 ghost cell' stencils)). */

    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    HERE is where we need to be VERY careful! As explained 
    [in the documentation](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README#a.-prediction-projection-steps) 
    for building a self-similar solver for the scale invariant problem of 
    [Keller \& Miksis (1983)](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c), 
    because the advection velocity IS NOT divergence 
    free, additional source terms appear compared to the reference paper! 
    Indeed, since $\textcolor{Orchid}{
    \left(\mathbf{\overline{\Lambda}} 
    \cdot \boldsymbol{\nabla} \right) \overline{\mathbf{u}} 
  } =
  \textcolor{Orchid}{
    \boldsymbol{\nabla} \cdot \left(\mathbf{\overline{\Lambda}} 
    \otimes \overline{\mathbf{u}} \right) 
  }
  \textcolor{green}{
    + \dfrac{2 N_d}{3} \, \overline{\mathbf{u}}
  }$, then:

    $$
    \partial_\tau \overline{u}_d^n 
    + \boldsymbol{\nabla} \cdot \left(\overline{u}_d^n 
    \textcolor{Orchid}{\mathbf{\overline{\Lambda}}^n} \right) 
    = \overline{g}_d^n
    \quad \Leftrightarrow \quad 
    \partial_\tau \overline{u}_d^n 
    =  
    - \textcolor{Orchid}{\mathbf{\overline{\Lambda}}^n} \cdot 
    \boldsymbol{\nabla} \overline{u}_d^n 
    + \overline{g}_d^n 
    \textcolor{green}{ + \dfrac{2 N_d}{3} \, \overline{u}_d^n}    
    $$

    and consequently we need to solve:

    $$
    \overline{u}_{f,p,d}^{n+1/2} 
    = \overline{u}_d^n 
    + \dfrac{\Delta}{2} \left[
      1 - \dfrac{\Delta \tau}{\Delta} \, \textcolor{Orchid}{\overline{\lambda}_d^n}
    \right] \partial_d \overline{u}_d^n  
    - \dfrac{\Delta \tau}{2} \textcolor{Orchid}{\overline{\lambda}_{\perp d}^n} \,
    \partial_{\perp d} \overline{u}_d^n 
    + \dfrac{\Delta \tau}{2} \left[
      \overline{g}_d^n 
      \textcolor{green}{ + \dfrac{2 N_d}{3} \, \overline{u}_d^n }
    \right]
    $$
    */
    double f2 = f[i] 
      + ( (src[] + src[-1]) + (2.*Nd/3.)*(f[] + f[-1]) )*dt/4. 
      + s*(1. - s*un)*g.x[i]*Delta/2.;

    /**
For tangential components... */

    #if dimension > 1
    if (fm.y[i] && fm.y[i,1]) {
      double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
      double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
      f2 -= dt*vn*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i] && fm.z[i,0,1]) {
      double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
      double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
      f2 -= dt*wn*fzz/(2.*Delta);
    }
    #endif

    flux.x[] = f2*uf.x[]; 
  }
}

/**
The function below uses the *selfsim_tracer_fluxes* function to integrate the
advection equation, using an explicit scheme with timestep *dt*, for
each tracer in the list. */

void selfsim_advection (scalar * tracers, face vector u, double dt, double Nd, 
		scalar * src = NULL)
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    selfsim_tracer_fluxes (f, u, flux, dt, Nd, source);
#if !EMBED
    foreach()
      foreach_dimension()
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}
