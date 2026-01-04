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


scalar f1[], f2[], f3[];

#ifndef SUBGRID
#define SUBGRID 1
#endif



#include "vof_with_bl_tracer_ARCELOR.h"
//#include "vof2.h"


scalar * interfaces = {f1,f2,f3};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;

face vector alphav[];
scalar rhov[];
face vector muv[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  mu = muv;
}

#ifndef rho
//We define rho based on the sum of product of fluid fraction and density
#define rho(f1,f2,f3) (clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) > 0 ? (clamp(f3,0.,1.)*rho3 + clamp(f2,0.,1.)*rho2 + clamp(f1,0.,1.)*rho1)/(clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) : 0.0
#endif
//We define also rhoM to avoid zero value for rho
#define rhom(f1,f2,f3) clamp(rho(f1,f2,f3), min(min(rho1, rho2), rho3), max(max(rho1, rho2), rho3))

#ifndef mu
//We define mu based on the sum of product of fluid fraction and viscosity
#define mu(f1,f2,f3) (clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) > 0 ? (clamp(f3,0.,1.)*mu3 + clamp(f2,0.,1.)*mu2 + clamp(f1,0.,1.)*mu1)/(clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) : 0.0
#endif
//We define also muM to avoid zero value for mu
#define mum(f1,f2,f3) clamp(mu(f1,f2,f3), min(min(mu1, mu2), mu3), max(max(mu1, mu2), mu3))

#ifdef FILTERED
scalar sf1[], sf2[], sf3[];
#else
# define sf1 f1
# define sf2 f2
# define sf3 f3
#endif



event tracer_advection(i++) {

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

  sf1.dirty = true;
  sf2.dirty = true;
  sf3.dirty = true;
#endif

}


event properties (i++)
{
  #if SUBGRID
    vector n1[], n2[], n3[];
    scalar alpha1[], alpha2[], alpha3[];

    if (f1.height.x.i){
      heights (f1, f1.height);
    }
    reconstruction (f1, n1, alpha1);

    if (f2.height.x.i){
      heights (f2, f2.height);
    }
    reconstruction (f2, n2, alpha2);

        if (f3.height.x.i){
      heights (f3, f3.height);
    }
    reconstruction (f3, n3, alpha3);
  #endif

  foreach_face() {

    #if SUBGRID
    double ff1 = 0;
    double ff2 = 0;
    double ff3 = 0;
    for (int i = -1; i <= 0; i++){
      ff1 += (f1[i] <= 0. || f1[i] >= 1.) ? f1[i]/2 :
      rectangle_fraction ((coord){2*((double)i + 0.5)*n1.x[i], n1.y[i], n1.z[i]}, alpha1[i],
        (coord){-0.5, -0.5, -0.5},
        (coord){0, 0.5, 0.5})/2;

      ff2 += (f2[i] <= 0. || f2[i] >= 1.) ? f2[i]/2 :
      rectangle_fraction ((coord){2*((double)i + 0.5)*n2.x[i], n2.y[i], n2.z[i]}, alpha2[i],
        (coord){-0.5, -0.5, -0.5},
        (coord){0, 0.5, 0.5})/2;

      ff3 += (f3[i] <= 0. || f3[i] >= 1.) ? f3[i]/2 :
      rectangle_fraction ((coord){2*((double)i + 0.5)*n3.x[i], n3.y[i], n3.z[i]}, alpha3[i],
        (coord){-0.5, -0.5, -0.5},
        (coord){0, 0.5, 0.5})/2;
    }
    #else
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    double ff3 = (sf3[] + sf3[-1])/2.;
    #endif

    alphav.x[] = (fm.x[]/rhom(ff1,ff2,ff3));
    muv.x[] = (fm.x[]*mum(ff1,ff2,ff3));
  }

  foreach(){
    rhov[] = (cm[]*rhom(sf1[],sf2[],sf3[]));
  }

#if TREE
  sf1.prolongation = fraction_refine;
  sf2.prolongation = fraction_refine;
  sf3.prolongation = fraction_refine;
  sf1.dirty = true;
  sf2.dirty = true;
  sf3.dirty = true;
#endif
}
