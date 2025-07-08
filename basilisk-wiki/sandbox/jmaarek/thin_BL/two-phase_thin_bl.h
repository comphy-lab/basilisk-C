/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */


/*We use a sharper calculation of density and viscosity at the interface. 
This modifies the local velocity field at the interface which is required for accurate calculation of mass transfer*/

#define SUBGRID 1
#include "vof_with_bl_tracer.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
#if SUBGRID
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
#endif
}

event properties (i++)
{
  #ifdef SUBGRID
    vector n[];
    scalar alpha[];

    if (f.height.x.i){
      heights (f, f.height);
    }
    reconstruction (f, n, alpha);
  #endif

  foreach_face(){
    #if SUBGRID
    double ff = 0;
    for (int i = -1; i <= 0; i++){
      ff += (f[i] <= 0. || f[i] >= 1.) ? f[i]/2 :
      rectangle_fraction ((coord){2*((double)i + 0.5)*n.x[i], n.y[i], n.z[i]}, alpha[i],
        (coord){-0.5, -0.5, -0.5},
        (coord){0, 0.5, 0.5})/2;
    }
    #else
    double ff = (sf[] + sf[-1])/2.;
    #endif
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff);
    }
  }
  
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE  
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
