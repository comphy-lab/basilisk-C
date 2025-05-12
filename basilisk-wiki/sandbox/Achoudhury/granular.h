#include "navier-stokes/centered.h"
#include "vof.h"
#define RHOF 1e-3
#define mug  1e-4

face vector alphav[];
face vector muv[];
face vector eta[];
scalar rhov[];
// note the v
scalar f[];
scalar * interfaces = {f};

double Co=0,D=1./30,mus=0.32,dmu=.28,I0=.4;
double etamax=100;

#define rho(f) ((f) + RHOF*(1. - (f)))

event init_granul (t = 0) {
    alpha = alphav;
    mu = muv;
    rho = rhov;
}

event properties (i++) {
    
    trash ({alphav});
    scalar eta[];
    
    scalar fa[];
    foreach()
    fa[] = (4.*f[] +
            2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
            f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
    
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
                double In = D2*D/sqrt(p[]);
                double muI = mus + dmu*In/(I0 + In);
                double etamin = sqrt(D*D*D);
                eta[] = max((Co + muI*p[])/D2, etamin);// this D2 is sqrt(2) D2
                eta[] = min(eta[],etamax);      }
        }
    }
    boundary ({eta});

    foreach_face() {
        double fm = (fa[] + fa[-1,0])/2.;
        muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
        // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug); was in Gerris
        alphav.x[] = 1./rho(fm);
    }
    foreach()
    rhov[] = rho(fa[]);
    boundary ({muv,alphav,rhov});
}
