#include "diffusion.h"
/** In progress. To be improved : correct interpolation of h_f and area_f.*/

/**# Tracer diffusion inside a layer*/

double D_hor = 0.;

/**
We consider the diffusion of a tracer $tr$ within each layer. As the thickness $h_l$ is not constant, it stays The diffusion equation within the layer is discretized :

$$
h_l \frac{\partial tr_l}{\partial t} = - \frac{\partial}{\partial x} J\, h_l \quad \text{with} 
\quad J = - D\, \nabla_l tr_l \quad
\text{and} 
\quad \nabla_l \approx \frac{\Delta}{\mathcal{A}}\, \frac{\partial }{\partial x}
$$
$$
\text{then, with a 1D gradient} \quad \, h_l\, \frac{\partial tr_l}{\partial t} =
\nabla(\frac{D\, h_l\,}{\mathcal{A}} \nabla tr_l)
$$
Therefore we see that to embody horizontal diffusion within the layer,
we need to :

* multiplie by $h_l$ the left term and the flux term, to take into account the fact that the thickness depends on x
* divid by the interface length the diffusion coefficient to take into acount the curvature of the layer. */

/**# Diffusion over the surface
In the same way, surfactant diffusion over the surface can be written as 
$$
\mathcal{A}\, \frac{\partial \Gamma}{\partial t} =
\nabla(\frac{D\,}{\mathcal{A}} \nabla \Gamma)
$$*/


/**# Diffusion function
The inputs of the function are:

* $tr$: diffusive tracer field,
* $dt$: the time step,
* $h$: layer thickness field,
* $area$: area ratio field (compared to $\Delta$),
* $surface$ boolean two select between surface or in-layer diffusion.
 */
  
mgstats horizontal_diffusion (scalar tr, double dt, scalar h, scalar area, bool surface)
{
  /**
  We allocate the face fields for the weighted diffusion coefficient and the theta correction. */
  face vector diffusion_coefficient[];
  scalar theta_c[], h_f[], area_f[];
  
  /**
  The theta field is computed : either with area or height field, depending on the boolean. Then the diffusion coefficient is computed.*/
  foreach()
    theta_c[] = surface ? area[]: h[];
  boundary({theta_c});
      
  foreach_face(){
    h_f[] = (h[]+h[-1])/2.;
    area_f[] = (area[]+area[-1])/2.;
    diffusion_coefficient.x[] = D_hor*h_f[]/area_f[];
  }
  boundary((scalar*) {diffusion_coefficient});
/**
  The implicit diffusion solver is called.*/
  return diffusion (tr, dt, D = diffusion_coefficient, theta = theta_c);
}


/**
# Horizontal diffusion event

The function, in the layers, is called during the event diffusive_term, just after the viscous_term event.*/

event diffusive_term (i++,last)
{  
  if (D_hor > 0.){
    scalar area_l[], zl[];
    foreach()
      zl[] = 0;
    scalar h, c;
    
   for (h, c in hl, cl){
      foreach()
        zl[] += h[]/2/Delta;
      boundary ({zl});

      foreach()
        area_l[] = sqrt(1. + sq((zl[1] - zl[-1])*0.5));
      
      foreach()
        zl[] += h[]/2/Delta;
        
      horizontal_diffusion (c, dt, h, area_l, false);
      boundary ({c});
    }
  }
}
