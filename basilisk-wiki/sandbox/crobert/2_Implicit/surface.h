/**
# Surface species
The file surface.h aims to add surface species in the multilayer solver.
These species will be advected on the surface using a BCG scheme and
 diffused on the surface if a diffusion coefficient is given as an
 attribute in init event.*/
double D_surf = 0.;  //Better if added as an attribute (for other fields)

/** 
To ensure mass conservation and independance toward the area metric $\mathcal{A}$, 
the surface concentration field must be :
$$
M = \Gamma * \mathcal{A}
$$
We define a function for area computation.
*/

double area (Point point) {
  double A = 1.;
  foreach_dimension ()
    A *= (sqrt(1. + sq((eta[1] - eta[-1])/(2.*Delta))));
  return A;
}


/**
By default, a unique surface concentration field $M$ is defined.
 Other surface fields can be added using surface list.*/
scalar * surface = NULL;
scalar M[];
event defaults (i = 0)
{
  surface = list_append (surface, M);
}

/**
# Advection over the surface 
*/


void advect_surface (scalar * surface, face vector hu, face vector hf, double dt)
{ 
    /**
    We compute the flux $(Mu)_{i+1/2,k}$ for each tracer $M$, using a
    variant of the BCG scheme. */
    
  face vector flux[];
  for (scalar M in surface) 
    foreach_face() {
    	double un = dt*hu.x[0,0,nl-1]/(hf.x[0,0,nl-1]*Delta + dry), a = sign(un);
    	int i = -(a + 1.)/2.;
      double h_metric = sqrt(1. + sq((eta[] - eta[-1])/Delta));
    	double g = M.gradient ?
  	  M.gradient (M[i-1], M[i], M[i+1])/Delta :
  	  (M[i+1] - M[i-1])/(2.*Delta*h_metric);
    	double s2 = M[i] + a*(1. - a*un)*g*Delta/2.;
  	  flux.x[] = s2*hu.x[0,0,nl-1]/(hf.x[0,0,nl-1] + dry);
    }
    boundary_flux ({flux});

    /**
    We compute $(M)^{n+1}_i = (M)^n_i + \Delta t 
    [(Mu)_{i+1/2} -(Mu)_{i-1/2}]/\Delta$. */
    
   foreach()
	  foreach_dimension()
	    M[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);

  scalar * list = list_copy (surface);
  boundary (list);
  free (list);
}

#if IMPLICIT_H
event half_advection (i++) {
  if (theta_H < 1.)
    advect_surface (surface, hu, hf, (1. - theta_H)*dt);
}
#endif

event remap (i++) {
  advect_surface (surface, hu, hf, dt);
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
    theta_c[] = area(point);
    tr_l[] = tr[];
  }
  boundary({theta_c, tr_l});

  /**
  The diffusion coefficient is computed on each face.*/
  foreach_face(){
    double h_metric = sqrt(1. + sq((eta[] - eta[-1])/Delta));
    diffusion_coefficient.x[] = D/h_metric;
    }
  boundary((scalar*) {diffusion_coefficient});
  
  /**
  The implicit diffusion solver is called.*/
  mgstats s = diffusion (tr_l, dt, D = diffusion_coefficient, theta = theta_c);
  foreach()
    tr[] = tr_l[];
  delete ({theta_c, tr_l, diffusion_coefficient});
  return s;
}

event viscous_term (i++) {

  for (scalar M in surface)
    if(D_surf > 0) {
      foreach()
        M[] *= (1./area(point));
      boundary ({M});
  
      surface_diffusion (h, M, D_surf, dt);
      
      foreach()
        M[] *= area(point);  
      boundary ({M});
    }
}

/**
# Cleanup*/
event cleanup (t = end)
{
  free (surface), surface = NULL;
}