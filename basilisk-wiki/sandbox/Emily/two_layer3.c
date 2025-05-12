/** 
#Two layer model - Do nothing test 3
Test case with parabolic basin. Ensures if initial condition is just bottom layer with flat water on top that nothing happens.
*/
#include "two_layer.h"
#define MAXLEVEL 8
#define TOL 1e-8
scalar hinit[],hsinit[];

int main() {
  N = 1 << MAXLEVEL;
  run();
}

event init (i=0){
  foreach(){
    zb[]=sq(x)+sq(y)-sq(1.);
    eta[]=0.;
    hs[]=max(0.,-0.5 - zb[]);
    hsinit[]=hs[];
    h[]=max(0.,eta[] - zb[] - hs[]);
    hinit[]=h[];
  }
  boundary({zb,eta,hs,h});
}

event end (i=1000){
  foreach(){
    hinit[]-=h[];
    hsinit[]-=hs[];
  }
  assert(normf(hinit).max<TOL);
  assert(normf(hsinit).max<TOL);
}
    
