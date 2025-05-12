/**
#All mach scheme
*/

#include "allmach_weno_poisson4_rk4.h"

scalar zb[],h[];
double G=1.;
double dry=1e-4;

scalar rhoc2v[];
face vector alphav[], av[];

event defaults (i=0) {
  rho = h;
  alpha = alphav;
  a = av;
  rhoc2 = rhoc2v;

  foreach()
    h[] = 1.;
  boundary({h});

}

event init (i=0) {
  foreach()
    p[] = G*sq(h[])/2.;
  boundary({zb,p});
}

double CFLa = HUGE;
 
event stability (i++){
  if(CFLa < HUGE)
     foreach()
        if(h[] > dry){
           double dt = CFLa*Delta/sqrt(G*h[]);
           if(dt<dtmax)
              dtmax=dt;
          }
}


void properties (scalar hp, scalar rhoc2p, scalar psp, face vector alphap, face vector mup){

 foreach(){
   rhoc2p[] = G*sq(max(hp[],dry));
   psp[] = rhoc2p[]/2.;
 }

 foreach_face(){
   if ((hp[]   > dry && hp[-1] > dry)||
       (hp[]   > dry && hp[] + zb[] >= zb[-1])||
       (hp[-1] > dry && hp[-1] + zb[-1] >= zb[]))
           alphap.x[] = 2./(max(h[],dry)+max(h[-1],dry));
   else
           alphap.x[] = 0.;
 }

 boundary ({hp,rhoc2p,psp});
 boundary ((scalar *){alphap,mup});
  
}

void acceleration (face vector alphap, face vector ap){

  foreach_face(){
    if(alphap.x[])
      ap.x[] = -G*(zb[]-zb[-1])/Delta;
    else
      ap.x[] = 0.;
  }
  boundary((scalar *){ap});
 
}
