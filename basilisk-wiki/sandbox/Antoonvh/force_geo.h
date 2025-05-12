/**
# Pressure-gradient force in a rotating frame of reference

This file implements both a Coriolis force ($F_{c} = \mathbf{f} \times
\mathbf{v}$) and and a large-scale pressure-gradient force $\nabla_h p$.
Assuming that $\mathbf{f}$ points in the vertial direction (zenith)
with $f_y$= `f_cor` $>0$. Then we can associate a so-called
geostrophic velocity vector with the aforementionted focings.

$$\mathbf{U_g} = \frac{e_z}{\rho f}\times\nabla_h p$$
 
and $\mathbf{U_g} = \{$`U_GEO`, `V_GEO`, 0$\}$, which are two macros
that should be `#define`d before including this file, which should
happen later than the inclusion of the [centered
slover](/src/navier-stokes/centered.h).

note that the meteorological convention is rotated: $\{u,v,w\} =$
`{u.z, u.x, u.y}. 
*/
#ifndef U_GEO
#define U_GEO 0.      
#endif
#ifndef V_GEO
#define V_GEO 0.
#endif
double  f_cor = 1e-4; //A non-zero constant
/**
## Implementation

in case of a two dimensional setup, we add a perpendicular velocity
field for the above equations to make sense. The positive direction is
$e_z = e_x \times e_y$ (towards you in typical visualizations).
*/
#if dimension == 2
#include "diffusion.h"
#include "tracer.h"
scalar uz[];
event tracer_diffusion (i++)
  diffusion (uz, dt, mu);
#endif

/**
   Inspired by other pages, we possibly need to allocate the
   acceleration field in an early event. 
 */
event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0;
    boundary ((scalar*){a});
#if dimension == 2
    tracers = list_append (tracers, uz);
#endif
  }
}
/**
   The "geostrophic velocity" gives a convinient basis to implement
   the pressure gradient force. 
*/
event acceleration (i++) {
  face vector av = a;
  boundary (all);
#if dimension == 3   // u.z should be available
  foreach_face(x) 
    av.x[] += -f_cor*((u.z[] + u.z[-1])/2.     - U_GEO);
  foreach_face(z) 
    av.z[] +=  f_cor*((u.x[] + u.x[0,0,-1])/2. - V_GEO);
#else                // we use uz
  foreach_face(x) 
    av.x[] +=   -f_cor*((uz[] + uz[-1])/2. - U_GEO);
  foreach()
    uz[]   += dt*f_cor*(u.x[]              - V_GEO);
#endif
}
