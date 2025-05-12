/**
# Three-phase interfacial flows

This file helps setup simulations for flows of three fluids separated by
corresponding interfaces (i.e. immiscible fluids). It can typically be used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h)
and [tension_three-phase.h](http://basilisk.fr/sandbox/vatsal/ThreePhase/tension_three-phase.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The method developed here is inspired by the discussions [here](https://groups.google.com/forum/#!topic/basilisk-fr/Sj_qNMXOfXE)

We use two different VOF tracers to track three fluids. The details of which is given in the following figure:
<figure>
<p align="center">
  <img src="https://www.dropbox.com/s/hiqxbu45hzhihm1/Schematic.png?dl=1" width="40%">
  <figcaption><p align="center">Figure 1. Schematic of the problem: Details of the VOF fields.</figcaption>
</figure>
So this method is similar to the one suggested by Jose. The difference is in the way we include the
surface tension force (described in [tension_three-phase.h](http://basilisk.fr/sandbox/vatsal/ThreePhase/tension_three-phase.h)).<br/>
The densities and dynamic viscosities for fluid 1, 2 and 3 are *rho1*,
*mu1*, *rho2*, *mu2*, and *rho3*, *mu3* respectively. */

#include "vof.h"
double VOF_cutoff = 0.01;
scalar f1[], f2[], *interfaces = {f1, f2};
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

  if (mu1 || mu2)
    mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic).<br/>
The properties (based on Figure 1) can be written as:
$$ \hat{A} (f_1, f_2) = f_1(1-f_2)\,A_l + f_1f_2\,A_2 + (1-f_1)\,A_3 \: \: \:  \forall  \: A \in \{\rho, \mu\} $$
The reason we choose this equation because the sum of coefficients of $A_i$ is zero.
*/

#ifndef rho
#define rho(f1, f2) (clamp(f1*(1-f2), 0., 1.) * rho1 + clamp(f1*f2, 0., 1.) * rho2 + clamp((1-f1), 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2) (clamp(f1*(1-f2), 0., 1.) * mu1 + clamp(f1*f2, 0., 1.) * mu2 + clamp((1-f1), 0., 1.) * mu3)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event properties (i++) {
  /**
  I do not understand this part of the code. In principle, f2 > 0.5 and f1 < 0.5
  should not occur because of the way they are described.
  But I see it happening sometimes (usually at triple contact point).
  Hence, I include this if-else loop.
  */
  if (i > 1){
    foreach(){
      if (f2[] > 0.5 & f1[] < 0.5){
        f1[] = f2[];
      }
    }
  }

  /**
  When using smearing of the density jump, we initialise sf_i with the
  vertex-average of f_i. */
#ifdef FILTERED
  int counter1 = 0;
  for (scalar sf in smearInterfaces){
    counter1++;
    int counter2 = 0;
    for (scalar f in interfaces){
      counter2++;
      if (counter1 == counter2){
        // fprintf(ferr, "%s %s\n", sf.name, f.name);
      #if dimension <= 2
          foreach(){
            sf[] = (4.*f[] +
        	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
        	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
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
          }
      #endif
      }
    }
  }
  #endif
#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = refine_bilinear;
    boundary ({sf});
  }
#endif

  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff1, ff2);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1, ff2);
  }
  foreach(){
    rhov[] = cm[]*rho(sf1[], sf2[]);
  }

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    boundary ({sf});
  }
#endif
}
