/**
# Surface species
The file surface.h aims to add surface species in the multilayer solver.
These species will be advected on the surface using a BCG scheme and diffused 
on the surface if the diffusion coefficient is given as attributes in init event.
 
To ensure mass conservation and independance toward the area metric $\mathcal{A}$, 
the surface concentration field must be :
$$
M = \Gamma * \mathcal{A}
$$
By default, a unique surface concentration field $M$ is defined. 
Other surface fields can be added using surface list and the diffusion coefficient attribute.*/
#include "curvature.h"
scalar * surface = NULL;

attribute {
  double D_surf;
}

scalar M[];
event defaults (i = 0)
{
  surface = list_append (surface, M);
}

/**
# Advection over the surface */
event advection_term (i++)
{ 
  for (scalar M in surface) {
    face vector flux[];
    foreach_face() {
      double un = 0;
      if (eta[0,0,nl-1] > dry)
        un = hu.x[0,0,nl-1]/hf.x[0,0,nl-1];
        un += (dut.x[] + dut.x[-1])/2.*hf.x[0,0,nl-1]/2.;
      double a = ((un) > 0 ? 1 : -1);
  	  int i = -(a + 1.)/2.;
      double g = (M[i+1] - M[i-1])/(2.*hm(point, false));
  	  double s2 = M[i] + a*(1. - a*dt*un/Delta)*g/2.;
      flux.x[] = s2*un;
    }
    boundary_flux ({flux});
  
    foreach()
      foreach_dimension ()
        M[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[0,0,nl-1]);
    boundary({M});
  }
}

/**
# Diffusion over the surface
Surfactant diffusion over the surface can be written as 
$$
\mathcal{A}\, \frac{\partial \Gamma}{\partial t} =
\nabla(\frac{D\,}{\mathcal{A}} \nabla \Gamma)
$$*/
#include "diffusion.h"
mgstats surface_diffusion (scalar h, scalar tr, double D, double dt)
{
  /**
  We allocate the tracer and parameters fields.*/
  vector diffusion_coefficient = new face vector[nl];
  scalar theta_c = new scalar[nl];
  scalar tr_l = new scalar[nl];
  reset({theta_c, tr_l, diffusion_coefficient}, 0);

  /**
  In each layer, the tracer is copied and the theta field is computed as the cell height.*/
  foreach(){
    theta_c[] = area_m(point);
    tr_l[] = tr[];
  }
  boundary({theta_c, tr_l});

  /**
  The diffusion coefficient is computed on each face.*/
  foreach_face()
    diffusion_coefficient.x[] = D/hm(point, false);
  boundary((scalar*) {diffusion_coefficient});
  
  /**
  The implicit diffusion solver is called.*/
  mgstats s = diffusion (tr_l, dt, D = diffusion_coefficient, theta = theta_c);
  foreach()
    tr[] = tr_l[];
  delete ({theta_c, tr_l, diffusion_coefficient});
  return s;
}

event viscous_term (i++)
{
  for (scalar M in surface)
    if(M.D_surf) {
      foreach()
        M[] *= (1./area_m(point));
      boundary ({M});
  
      surface_diffusion (h, M, M.D_surf, dt);
      
      foreach()
        M[] *= area_m(point);  
      boundary ({M});
    }  
}

event cleanup (t = end)
{
  free (surface), surface = NULL;
}