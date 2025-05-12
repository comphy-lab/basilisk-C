/**
# Three-phase interfacial flows

This file helps setup simulations for flows of three fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid method.   
The volume fraction in fluid 2 is **f2=1**, in fluid 3 is **f3=1**, and in fluid 1 is **f1=(1 -f2 -f3)**.   
The densities and dynamic viscosities for fluid 1, 2, and 3 
are *rho1*, *mu1*, *rho2*, *mu2*, *rho3*, *mu3*, respectively. */

#include "vof.h"

scalar f2[],f3[],f1[];
scalar * interfaces = {f2,f3,f1};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;

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
  
  if (mu1 || mu2 || mu3)
    mu = new face vector;
}

/**
Now we define new relations for rho et mu with fields f2, f3, and f1.   
note that we must have f1 + f2 + f3 = 1.0 */

#ifndef rho
// two phases definition was rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
// three phases definition
#define rho(f2,f3) (clamp(1.-f2-f3,0,1)*rho1+clamp(f2,0,1)*rho2+clamp(f3,0,1)*rho3)
#endif

#ifndef mu
// two phases definition was mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
// three phases definition
#define mu(f2,f3)  (clamp(1.-f2-f3,0,1)*mu1+clamp(f2,0,1)*mu2+clamp(f3,0,1)*mu3)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf1[],sf2[],sf3[];
#else
# define sf1 f1
# define sf2 f2
# define sf3 f3
#endif

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf1
#if dimension <= 2
  foreach()
    sf1[] = (4.*f1[] + 
	    2.*(f1[0,1] + f1[0,-1] + f1[1,0] + f1[-1,0]) +
	    f1[-1,-1] + f1[1,-1] + f1[1,1] + f1[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf1[] = (8.*f1[] +
	    4.*(f1[-1] + f1[1] + f1[0,1] + f1[0,-1] + f1[0,0,1] + f1[0,0,-1]) +
	    2.*(f1[-1,1] + f1[-1,0,1] + f1[-1,0,-1] + f1[-1,-1] + 
		f1[0,1,1] + f1[0,1,-1] + f1[0,-1,1] + f1[0,-1,-1] +
		f1[1,1] + f1[1,0,1] + f1[1,-1] + f1[1,0,-1]) +
	    f1[1,-1,1] + f1[-1,1,1] + f1[-1,1,-1] + f1[1,1,1] +
	    f1[1,1,-1] + f1[-1,-1,-1] + f1[1,-1,-1] + f1[-1,-1,1])/64.;
#endif
#endif

#ifndef sf2
#if dimension <= 2
  foreach()
    sf2[] = (4.*f2[] + 
	    2.*(f2[0,1] + f2[0,-1] + f2[1,0] + f2[-1,0]) +
	    f2[-1,-1] + f2[1,-1] + f2[1,1] + f2[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf2[] = (8.*f2[] +
	    4.*(f2[-1] + f2[1] + f2[0,1] + f2[0,-1] + f2[0,0,1] + f2[0,0,-1]) +
	    2.*(f2[-1,1] + f2[-1,0,1] + f2[-1,0,-1] + f2[-1,-1] + 
		f2[0,1,1] + f2[0,1,-1] + f2[0,-1,1] + f2[0,-1,-1] +
		f2[1,1] + f2[1,0,1] + f2[1,-1] + f2[1,0,-1]) +
	    f2[1,-1,1] + f2[-1,1,1] + f2[-1,1,-1] + f2[1,1,1] +
	    f2[1,1,-1] + f2[-1,-1,-1] + f2[1,-1,-1] + f2[-1,-1,1])/64.;
#endif
#endif 

#ifndef sf3
#if dimension <= 2
  foreach()
    sf3[] = (4.*f3[] + 
	    2.*(f3[0,1] + f3[0,-1] + f3[1,0] + f3[-1,0]) +
	    f3[-1,-1] + f3[1,-1] + f3[1,1] + f3[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf3[] = (8.*f3[] +
	    4.*(f3[-1] + f3[1] + f3[0,1] + f3[0,-1] + f3[0,0,1] + f3[0,0,-1]) +
	    2.*(f3[-1,1] + f3[-1,0,1] + f3[-1,0,-1] + f3[-1,-1] + 
		f3[0,1,1] + f3[0,1,-1] + f3[0,-1,1] + f3[0,-1,-1] +
		f3[1,1] + f3[1,0,1] + f3[1,-1] + f3[1,0,-1]) +
	    f3[1,-1,1] + f3[-1,1,1] + f3[-1,1,-1] + f3[1,1,1] +
	    f3[1,1,-1] + f3[-1,-1,-1] + f3[1,-1,-1] + f3[-1,-1,1])/64.;
#endif
#endif

#if TREE
  sf1.prolongation = refine_bilinear;
  sf2.prolongation = refine_bilinear;
  sf3.prolongation = refine_bilinear;    
  boundary ({sf1,sf2,sf3});
#endif
  
  foreach_face() {
//    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    double ff3 = (sf3[] + sf3[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff2,ff3);
    if (mu1 || mu2 || mu3) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff2,ff3);
    }
  }
  foreach()
    rhov[] = cm[]*rho(sf2[],sf3[]);

#if TREE  
  sf1.prolongation = fraction_refine;
  sf2.prolongation = fraction_refine;
  sf3.prolongation = fraction_refine;
  boundary ({sf1,sf2,sf3});
#endif
}