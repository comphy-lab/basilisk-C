/**
# Centered phase change
This file extends the centered Navier-Stokes solver to account for phase change
in two-phase flows with porous media.
It modifies the projection method to include a gas source term in the continuity equation,
and it provides an advection scheme that can account for porosity variations.
*/

#include "poisson.h"
extern double rhoG;
extern scalar porosity, f;

scalar gas_source[];
scalar drhodt[];

/**
## Projection method with gas source term
We modify the Projection method to account for the gas source term
in the continuity equation.
*/

trace
mgstats project_sf (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;

    /**
    We add the gas source term to the divergence field
    */

    div[] += gas_source[]/dt;
#ifndef NO_EXPANSION
    div[] += drhodt[]/dt;
#endif
  }

#ifdef POROUS_ADVECTION
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);

  face vector alpha_eff[];
  foreach_face()
    alpha_eff.x[] = face_value(eps, 0)*alpha.x[];
#else
  #define alpha_eff alpha
#endif

  mgstats mgp = poisson (p, div, alpha_eff,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    uf.x[] -= dt*alpha_eff.x[]*face_gradient_x (p, 0);

  return mgp;
}

/**
## Advection with non diverging velocity field
We provide an advection scheme that can account for a diverging velocity field.
*/

#include "utils.h"
#include "bcg.h"

void advection_div (scalar * tracers, face vector u, double dt,
    scalar * src = NULL)
{
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
    tracer_fluxes (f, u, flux, dt, source);
#if !EMBED
    foreach() {
#if NO_ADVECTION_DIV
      double fold = f[];
#endif
      foreach_dimension()
#if NO_ADVECTION_DIV
        f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
#else
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#endif
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}


/**
## Default events and overrides
We set default values for the gas source and drhodt fields
*/

event defaults (i = 0) {
  foreach(){
    gas_source[] = 0.;
    drhodt[] = 0.;
  }
}

/** 
We set placeholder to set the correct order of the events 
*/

event set_dtmax (i++, last);
event stability (i++, last);
event reset_sources (i++, last);
event chemistry (i++, last);
event phasechange (i++, last);

/** 
We overwrite the project and advection events with the one defined above
*/
#define project(...) project_sf(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered.h"
#undef advection
#undef project


/**
# Porous media advection
We have the option to account for porous media advection by taking into 
account the porosity field in the advection term.
We set stokes=true to suppress the original advection term performed in the
centered.h file.
*/

#ifdef POROUS_ADVECTION
event defaults (i = 0) {
  stokes = true;
}

event advection_term (i++, last) {
  prediction();
  mgpf = project_sf (uf, pf, alpha, dt/2., mgpf.nrelax);

  /**
  porosity is a tracer field appended to f. Here we need the one-field form
  computed in the field 'eps'.
  */

  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);

  face vector ufn[];
  foreach_face() {
    double ef = face_value(eps, 0);
    ufn.x[] = uf.x[]/ef;
  }
  
  advection ((scalar *){u}, ufn, dt, (scalar *){g});
}

/** 
the stability event gets disable if stokes is set to true
since the solution is implicit. Therefore, we redefine the
stability event to set the timestep.
*/

event stability (i++, last) {
  dt = dtnext (timestep (uf, dtmax));
}
#endif
