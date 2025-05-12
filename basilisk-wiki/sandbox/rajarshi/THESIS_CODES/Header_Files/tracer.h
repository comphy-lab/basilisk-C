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
  #if WENO
     #if LIMITED
        #include "Adapt_Limiters.h"
     #else 
        #include "Adapt_No_Limiters.h"
     #endif  
     event defaults (i = 0) {
         for (scalar s in tracers) {
             s.refine = s.prolongation = refine_order5;
             s.restriction = restriction_volume_average;
         } 
     }
  #else // !WENO
     event defaults (i = 0) {
        for (scalar s in tracers) {
            s.refine  = refine_linear;
            s.restriction = restriction_volume_average;
        }
     }
  #endif // !WENO
#endif // !TREE

/**
The integration is performed using the Bell--Collela--Glaz scheme or the
5th-order WENO scheme. */

#if WENO
   #if LIMITED
      #include "WENO_Flux.h"
   #else
      #include "O5_Flux.h"
   #endif

   #include "runge-kutta.h"
   event tracer_advection (i++,last) {
        #if ORDER
          runge_kutta_D (tracers, t , dt, tracer_fluxes, ORDER);
        #else
          runge_kutta (tracers, t, dt, tracer_fluxes, 4);
        #endif
   }

#else // !WENO

    #include "bcg.h"
    event tracer_advection (i++,last) {
         advection (tracers, uf, dt);
    }

#endif

/**
Diffusion can be added by overloading this hook. */

event tracer_diffusion (i++,last);
