/**
# Tracer advection event

This event integrates advection equations of the form
$$
\partial_tf_i+\mathbf{u_f}\cdot\nabla f_i=0
$$
where $\mathbf{u_f}$ is the velocity field and $f_i$ are a list of
passive tracers.

The `tracers` list is defined elsewhere (typically by the user), the
face vector field `uf` and the timestep `dt` are defined by a
solver. */

extern scalar * tracers;
extern face vector uf;
extern double dt;

/**
On adaptive meshes, tracers need to use linear interpolation (rather
than the default bilinear interpolation) to ensure conservation when
refining cells. */

#if TREE
event defaults (i = 0) {
  for (scalar s in tracers) {
#if EMBED
    s.refine = refine_embed_linear;
    set_prolongation (s, refine_embed_linear);
#else
    s.refine  = refine_linear;
#endif
    set_restriction (s, restriction_volume_average);
  }
}
#endif

/**
The integration is performed using the Bell-Collela-Glaz scheme. */

#include "bcg.h"

event tracer_advection (i++,last) {
  advection (tracers, uf, dt);
}

/**
VOF advection happens here. */

event vof (i++, last);

/**
Diffusion can be added by overloading this hook. */

event tracer_diffusion (i++,last);
