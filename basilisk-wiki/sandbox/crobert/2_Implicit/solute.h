/** 
# Solutes
The file solute.h aims to add a solute in the multilayer solver. The solute
 will be advected as any tracers (hydro.h) and diffused vertically with the
 following diffusion coefficient.*/

double D_vert = 0.; //Better if added as an attribute (for other fields)
double D_hor = 0.;

/** By default, a unique solute concentration field $c$ is defined. 
Other solutal fields can be added using tracers and the diffusion 
coefficient attributes.*/
scalar * solute = NULL;
scalar c;
event defaults (i = 0)
{
  assert (nl > 0);
  c = new scalar[nl];
  reset({c}, 0);

  if (!linearised)
    tracers = list_append (tracers, c);
  solute = list_append (solute, c);
}

/**
# Diffusion
Vertical diffusion is called during the viscous event. A boundary condition can 
 be added at the bottom.*/
(const) scalar c_b = zeroc;
//(const) scalar c_lbda_b = {HUGE};
event viscous_term (i++)
{
  for (scalar c in solute) {
    if (D_hor > 0)
      horizontal_diffusion ({c}, D_hor, dt); //Added here but should respect timestep limitation (see diffusion.h)
    if (D_vert > 0)
      foreach()
        foreach_dimension()
  	      vertical_diffusion(point, h, c, dt, D_vert, 0, c_b[], HUGE);
    boundary ({c});
  }
}

/**
# Evaporation
Evaporation is controlled by the velocity $v_e$. 
It must be initialized to add evaporation. 

To ensure stability, the timestep must be controlled.
 The event set_dtmax is called to ensure that:
$$
v_e\,dt < 0.3\,h_{nl-1}
$$ 
*/
(const) scalar v_e = zeroc;

event stability (i++) 
{
  double dt_max = HUGE;
  foreach()
    if (v_e[] != 0)
      dt_max = min(dt_max, 0.3*h[0,0,nl-1]/v_e[]);
  dtmax = min(dtmax, dt_max);
}

/**This function calculates the evaporation term.
$$
\frac{dh}{dt} = v_e\,\mathcal{A}
$$
As the solute does not evaporate, the concentration field must be corrected
 such as: 
$$
c_{old}\, h_{old} = c\, h
$$
Evaporation is called just before advection (and eta computation). */

event viscous_term (i++)
{
  foreach()
    if (v_e[] != 0) {
      double area = 1.;
      foreach_dimension() {
        double dh = (eta[1] - eta[-1])/(2.*Delta);
        area *= sqrt(1 + sq(dh));
      }
      double h_evap = v_e[] * area * dt;
      h[0,0,nl-1] -= h_evap;
      for (scalar c in solute)
        c[0,0,nl-1] *= (1. + h_evap/h[0,0,nl-1]);
    }
}

/**
# Cleanup*/
event cleanup (t = end)
{
  delete ({c});
  free (solute), solute = NULL;
}