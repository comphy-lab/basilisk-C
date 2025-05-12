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

extern scalar cs, csm1;
extern bool emerged;

/**
On adaptive meshes, tracers need to use linear interpolation (rather
than the default bilinear interpolation) to ensure conservation when
refining cells. */

#if TREE
event defaults (i = 0) {
  for (scalar s in tracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
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
Diffusion can be added by overloading this hook. */

event tracer_diffusion (i++,last);

/**
We update the tracer in the solid and emerged cells, using the
boundary condition at the embedded boundary. */

event advection_term (i++,last)
{
  /**
  We first make sure not to use any values in newly emerged cells by
  setting the flag *emerged* to false. */

  emerged = false;
  boundary (tracers);
  
  /**
  In the solid cells, we set the scalar field *s* to 0. */

  for (scalar s in tracers)
    foreach()
      if (cs[] <= 0.) 
	s[] = 0.;	
  
  /**
  In the emerged cells, we use the solid boundary condition to
  update the scalar field. */

  for (scalar s in tracers)
    foreach() {
      if (csm1[] <= 0. && cs[] > 0.) {

	// Normal emerged cell
	if (cs[] < 1.) {

	  // Cell centroid, barycenter and normal of the embedded fragment
	  coord c, b, n;
	  embed_geometry (point, &b, &n);
	  double alpha = plane_alpha (cs[], n);
	  plane_center (n, alpha, cs[], &c); // Different from line_area_center

	  // Boundary condition on the embedded boundary
	  bool dirichlet = true;
	  double sb = (s.boundary[embed] (point, point, s, &dirichlet));
	  if (!dirichlet) {
	    double coef = 0.;
	    sb = neumann_scalar (point, s, cs, n, b, sb, &coef);
	    // s[] is undefined here so we use the average over the neighbors
	    int navg = 0;
	    double savg = 0.;
	    foreach_neighbor(1)
	      if (cs[] > 0. && (emerged || csm1[] > 0.)) {
		navg += 1;
		savg += s[];
	      }
	    sb += coef*(savg/(navg + SEPS));
	  }

	  // Emerged cell value
	  s[] = embed_extrapolate (point, s, cs, n, c, sb);
	}

	// Pathological emerged cell (cs = 1)
	else {
	  int navg = 0;
	  double savg = 0.;
	  foreach_neighbor(1)
	    if (cs[] > 0. && (emerged || csm1[] > 0.)) {
	      navg += 1;
	      savg += s[];
	    }
	  s[] = savg/(navg + SEPS);
	}
      }
    }
  
  /**
  Before using the *boundary* function, we set the *emerged* flag to
  true to indicate that all emerged cells have been updated and can
  now be used. */
  
  emerged = true;
  boundary (tracers);
}
