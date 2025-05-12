/**
# Integral formulation for surface tension

This page has moved [here](/src/integral.h).

The surface tension $\sigma$ will be associated to each levelset
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */

attribute {
  double sigma;
}

event defaults (i = 0) {
  
  /**
  Surface tension is a source term in the right-hand-side of the
  evolution equation for the velocity of the [centered Navier--Stokes
  solver](navier-stokes/centered.h) i.e. it is an acceleration. If
  necessary, we allocate a new vector field to store it. */

  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  We then consider each VOF interface with an associated value of
  $\sigma$ different from zero and set the maximum timestep. */
  
  for (scalar c in tracers)
    if (c.sigma) {
      double dt = sqrt (rhom*cube(dmin)/(pi*c.sigma));
      if (dt < dtmax)
	dtmax = dt;
    }
}

/**
## Curvature

This function computes the curvature from the levelset function *d* using:
$$
\kappa = \frac{d^2_xd_{yy} - 2d_xd_yd_{xy} + d^2_yd_{xx}}{|\nabla d|^3}
$$
*/

static inline double curvature (Point point, scalar d) {
  double dx = (d[1] - d[-1])/2.;
  double dy = (d[0,1] - d[0,-1])/2.;
  double dxx = d[1] - 2.*d[] + d[-1];
  double dyy = d[0,1] - 2.*d[] + d[0,-1];
  double dxy = (d[1,1] - d[-1,1] - d[1,-1] + d[-1,-1])/4.;
  double dn = sqrt(sq(dx) + sq(dy)) + 1e-30;
  return (sq(dx)*dyy - 2.*dx*dy*dxy + sq(dy)*dxx)/cube(dn)/Delta;
}

#define LINEAR_CURVATURE 0 // set to 1 to use linear interpolation of curvature

/**
## Surface tension term

The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  
  /**
  We check whether surface tension is associated with any levelset
  function *d*. */

  for (scalar d in tracers)
    if (d.sigma) {

#if LINEAR_CURVATURE
      /**
      We first compute the curvature. */

      scalar kappa[];
      foreach()
	kappa[] = curvature (point, d);
      boundary ({kappa});
#endif
      
      /**
      We allocate the surface tension stress tensor. */
      
      tensor S[];

      /**
      We compute the diagonal components of the tensor. */
      
      foreach()
	foreach_dimension() {
	  S.y.y[] = 0.;
	  for (int i = -1; i <= 1; i += 2)
	    if (d[]*(d[] + d[i]) < 0.) {
	      double xi = d[]/(d[] - d[i]);
	      double nx = ((d[1] - d[-1])/2. +
			   xi*i*(d[-1] - 2.*d[] + d[1]))/Delta;
#if LINEAR_CURVATURE // does not make much difference
              double ki = kappa[] + xi*(kappa[i] - kappa[]);
#else
	      double ki = curvature (point, d);
#endif
	      S.y.y[] += d.sigma*(fabs(nx)/Delta - sign(d[])*ki*(0.5 - xi));
	    }
        }
      boundary ({S.x.x, S.y.y});

      /**
      We compute the off-diagonal components of the tensor.  */
      
      foreach_vertex()
	foreach_dimension() {
	  if ((d[] + d[0,-1])*(d[-1] + d[-1,-1]) > 0.)
	    S.x.y[] = 0.;
	  else {
	    double xi = (d[-1] + d[-1,-1])/(d[-1] + d[-1,-1] - d[] - d[0,-1]);
	    double ny = (xi*(d[] - d[-1] + d[-1,-1] - d[0,-1]) +
			 d[-1] - d[-1,-1])/Delta;
	    S.x.y[] = - d.sigma*sign(d[] + d[0,-1])*ny/Delta;
	  }
        }

      /**
      Finally, we add the divergence of the surface tension tensor to
      the acceleration. */
      
      face vector av = a;
      foreach_face()
	av.x[] += (S.x.x[] - S.x.x[-1] + S.x.y[0,1] - S.x.y[])/Delta;      
    }
}
