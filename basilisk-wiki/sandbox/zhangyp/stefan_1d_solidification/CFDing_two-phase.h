/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 
 */


#include "fractions.h"//This code only called fraction_refine defined in fractions.h

scalar * interfaces = NULL;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ and average viscosity $\mu$ (on faces) as well
as the cell-centered density. */

face vector alphav[], muv[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  mu = muv;
}

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */
  //double c1, c2, c3;
  foreach_face() {
    double phase1 = (phi1_ns[] + phi1_ns[-1] )/2.;
    double phase2 = 1.0-(phi1_ns[] + phi1_ns[-1] )/2.;

    phase1 = 	clamp (phase1, 0., 1.);
    phase2 = 	clamp (phase2, 0., 1.);

    alphav.x[] = fm.x[]/(phase1 + phase2*rd);
    muv.x[] = fm.x[]*(phase1 + phase2*rv)/Re;
  }
  foreach()
  { 
    double phase1 = phi1_ns[];
    double phase2 = 1.0-phi1_ns[];

    phase1 = 	clamp (phase1, 0., 1.);
    phase2 = 	clamp (phase2, 0., 1.);

    rhov[] = cm[]*( phase1 + phase2*rd);
  }
  #if TREE  
  phi1.refine=phi1.prolongation = fraction_refine;
  boundary ({phi1});
  #endif 

}