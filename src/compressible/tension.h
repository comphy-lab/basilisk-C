/**
# Surface tension effects for the compressible solver 

This module incorporates surface tension effects in the compressible
solver. Unlike momentum, which is defined for the averaged mixture,
the energy is defined separately for both phases. For that reason, we
need to reconstruct the pressure of each phase from the averaged
pressure obtained from the Helmholtz--Poisson equation.

To reconstruct $p_1$ and $p_2$ from the averaged pressure $p$,
we interpret the averaged pressure as
$$p = f p_1 + (1-f) p_2$$
Using the Laplace equation
$$p_2 - p_1 = -\sigma \kappa + [[\mu n \cdot \tau \cdot n]]$$
we can solve the system of the two equations for $p_1$ and $p_2$
$$  p_1  = p + (1 - f) \sigma \kappa - (1-f) [[\mu n \cdot \tau \cdot n]]$$
$$  p_2  = p - f \sigma \kappa + f [[\mu n \cdot \tau \cdot n]]$$
The flux terms depending on pressure entering into the conservative total energy equation are
$$\nabla \cdot (u p_1) = \nabla \cdot (u p) + \nabla \cdot (u (1-f) \sigma \kappa) - 
                         \nabla \cdot (u (1-f) [[\mu n \cdot \tau \cdot n]])$$
$$\nabla \cdot (u p_2) = \nabla \cdot (u p) - \nabla \cdot (u f \sigma \kappa) + 
                         \nabla \cdot (u f [[\mu n \cdot \tau \cdot n]])$$

In this module we only introduce the correction due to surface tension. */

event end_timestep (i++)
{
  if (f.sigma > 0.) {
    scalar skappa[];
    curvature (f, skappa, f.sigma);

    /** 
    Here we just introduce the correction due to the surface tension
    contribution to the pressure jump. */
    
    face vector upf[];
    for (scalar fE in {fE1, fE2}) {
      foreach_face() {
	double fr = f[1], fl = f[];
	if (!fE.inverse)
	  fr = 1. - fr, fl = 1. - fl;
	if (fr + fl > 0. && fr + fl < 2.) {
	  if (fr > fl)
	    upf.x[] = uf.x[]*fl*skappa[];
	  else
	    upf.x[] = uf.x[]*fr*skappa[1];
	}
	else
	  upf.x[] = 0.;
      }

      foreach() {
	double div = 0.;
	foreach_dimension()
	  div += upf.x[1] - upf.x[];
	fE[] -= (fE.inverse ? 1. - f[] : f[])*div/Delta*dt/cm[];
      }
    }
  }
}
