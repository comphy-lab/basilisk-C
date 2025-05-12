/**
# Granular rheology, $\mu(I)$

## Purpose
This is the implementation of $\mu(I)$ rheology for bulk flows (without vof)

## Implementation
 We use
*/
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
//#include "navier-stokes/perfs.h"

#define mug  1e-4

double Co=0,D=1./30,mus=0.32,dmu=.28,I0=.4;
double etamax=10000;


face vector muv[];
scalar In[];
scalar muI[];

event init_granul (t = 0) {
    mu = muv;
}

/**
## computing the viscosity
 */
event properties (i++) {

scalar eta[];

foreach() {
        eta[] = mug;
        if (p[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            if (D2 > 0.) {
                D2 = sqrt(2.*D2)/(2.*Delta); // this D2 is sqrt(2) D2
                //double In = D2*D/sqrt(p[]);
                //double muI = mus + dmu*In/(I0 + In);
                In[] = D2*D/sqrt(p[]);
                muI[] = (mus + dmu*In[]/(I0 + In[]));                
                double etamin = sqrt(D*D*D);
                eta[] = max((Co + muI[]*p[])/D2, etamin);// this D2 is sqrt(2) D2
                eta[] = min(eta[],etamax);     
            }
        }
    }
  boundary ({eta});
  foreach_face()
      //muv.x[] = fm.x[]*(eta[] + eta[-1,0])/2.;
      muv.x[] = fm.x[]*(4.*eta[] +
            2.*(eta[-1,0] + eta[1,0] + eta[0,-1] + eta[0,1]) +
            eta[1,1] + eta[-1,1] + eta[1,-1] + eta[-1,-1])/16.;
      //muv.x[] = fm.x[]*1.; // For Newtonian limit
  boundary ((scalar*){muv});
}