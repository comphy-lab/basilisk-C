/** 
#Two layer model Do Nothing Test 1
Test case with parabolic basin. Ensures if initial condition is just flat water with no bottom layer that nothing happens.
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
    hs[]=0.;
    h[]=max(0.,eta[] - zb[] - hs[]);
    hinit[]=h[];
  }
  boundary({zb,eta,hs,h});
}

event end (i=1000){
  assert(normf(hs).max<TOL);
  foreach()
    hinit[]-=h[];
  assert(normf(hinit).max<TOL);
}
