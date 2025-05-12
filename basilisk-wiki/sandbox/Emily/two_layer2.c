/** 
Two layer model - Do nothing test 2
Test case with parabolic basin. Ensures if initial condition is just flat bottom layer with no water that nothing happens.
*/
#include "two_layer.h"
#define MAXLEVEL 8
#define TOL 1e-8
scalar hinit[];

int main() {
  N = 1 << MAXLEVEL;
  run();
}

event init (i=0){
  foreach(){
    zb[]=sq(x)+sq(y)-sq(1.);
    eta[]=0.;
    h[]=0.;
    hs[]=max(0.,eta[] - zb[] - h[]);
    hinit[]=hs[];
  }
  boundary({zb,eta,hs,h});
}

event end (i=1000){
  assert(normf(h).max<TOL);
  foreach()
    hinit[]-=hs[];
  assert(normf(hinit).max<TOL);
}
    
