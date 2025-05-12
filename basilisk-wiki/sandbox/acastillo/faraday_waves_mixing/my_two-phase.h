/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/** 
***We treat fluid 1 as a binary mixture characterised by the normalized
concentration field $c \in [0,1]$, such that
$$ \rho_1 = \rho_{1,max} c + \rho_{1,min} (1-c) $$
which is fine as long as $c$ does not seep into fluid 2.***
*/
extern scalar c;
double rho1_max = 1., rho1_min = 1.;
#define rho1(c) ( rho1_min + clamp(c,0.,1.)*(rho1_max - rho1_min) )

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
# define rho(f,c) (clamp(f,0.,1.)*(rho1(c) - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
scalar cf[];
#else
# define sf f
# define cf c
#endif

event tracer_advection (i++)
{

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. ***We also initialise cs with the vertex-average of
  c *** 
  */

#ifndef sf
#if dimension <= 2
  foreach(){
    sf[] = (4.*f[] +
      2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
      f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
    cf[] = (4.*c[] +
      2.*(c[0,1] + c[0,-1] + c[1,0] + c[-1,0]) +
      c[-1,-1] + c[1,-1] + c[1,1] + c[-1,1])/16.;
  }
#else // dimension == 3
  foreach(){
    sf[] = (8.*f[] +
      4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
      2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
    f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
    f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
      f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
      f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
    cf[] = (8.*c[] +
      4.*(c[-1] + c[1] + c[0,1] + c[0,-1] + c[0,0,1] + c[0,0,-1]) +
      2.*(c[-1,1] + c[-1,0,1] + c[-1,0,-1] + c[-1,-1] +
    c[0,1,1] + c[0,1,-1] + c[0,-1,1] + c[0,-1,-1] +
    c[1,1] + c[1,0,1] + c[1,-1] + c[1,0,-1]) +
      c[1,-1,1] + c[-1,1,1] + c[-1,1,-1] + c[1,1,1] +
      c[1,1,-1] + c[-1,-1,-1] + c[1,-1,-1] + c[-1,-1,1])/64.;
  }
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
  cf.prolongation = refine_bilinear;
  cf.dirty = true; // boundary conditions need to be updated
#endif
}

/**
***Finally, we include the dependency on f and c to update $\rho$ and $1/\rho$*** 
*/

event properties (i++)
{
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    double cc = (cf[] + cf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff,cc);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff);
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[],cf[]);

#if TREE
  sf.prolongation = fraction_refine;
  cf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
  cf.dirty = true; // boundary conditions need to be updated
#endif
}

#undef rho1
