/**
# Three phase interfacial flows

This is a modified version of [two--phase](/src/two-phase.h)
which allow to add a third phase explicitly.

This file helps setup simulations for flows of three fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](/src/navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f1=1$, $f2=1$ in fluid
2 and $f3=1$ in fluid 3. The densities and dynamic viscosities for fluid 1,2 are *rho1*,
 *mu1*, *rho2*, *mu2*, *rho3*, *mu3* respectively. */

#include "vof.h"

scalar f1[], f2[], f3[];
scalar * interfaces = {f1,f2,f3};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;

face vector alphav[];
scalar rhov[];


event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  if (mu1 || mu2|| mu3)
    mu = new face vector;
}

#ifndef rho
//We define rho based on the sum of product of fluid fraction and density
#define rho(f1,f2,f3) (clamp(f3,0.,1.)*rho3 + clamp(f2,0.,1.)*rho2 + clamp(f1,0.,1.)*rho1) 
//We define also rhoM to avoid zero value for rho
#define rhom(f1,f2,f3) max(rho(f1,f2,f3),rho1)
#endif

#ifndef mu
//We define mu based on the sum of product of fluid fraction and viscosity
#define mu(f1,f2,f3) (clamp(f3,0.,1.)*mu3 + clamp(f2,0.,1.)*mu2 + clamp(f1,0.,1.)*mu1) 
//We define also muM to avoid zero value for mu
#define mum(f1,f2,f3) max(mu(f1,f2,f3),mu1)
#endif

#ifdef FILTERED
scalar sf1[], sf2[], sf3[];
#else
# define sf1 f1
# define sf2 f2
# define sf3 f3
#endif



event properties(i++) {

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

  boundary ({sf1, sf2, sf3});
#endif

  foreach_face() {
    double f1m = (sf1[] + sf1[-1])/2.;
    double f2m = (sf2[] + sf2[-1])/2.;
    double f3m = (sf3[] + sf3[-1])/2.;
    
    alphav.x[] = (fm.x[]/rhom(f1m,f2m,f3m));
    if (mu1 || mu2 || mu3) {
      face vector muv = mu;
      muv.x[] = (fm.x[]*mum(f1m,f2m,f3m));
    }
  }
  foreach(){
    rhov[] = (cm[]*rhom(sf1[],sf2[],sf3[]));
  }

#if TREE  
  sf1.prolongation = fraction_refine;
  sf2.prolongation = fraction_refine;
  sf3.prolongation = fraction_refine;
  boundary ({sf1, sf2, sf3});
#endif
}

