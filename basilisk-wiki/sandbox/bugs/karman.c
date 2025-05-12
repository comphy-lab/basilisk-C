/**
# Source file crashes the compiler

This is due to the formatting of vectors without spaces, see below. */

#include "navier-stokes/mac.h"
#include "tracer.h"

/*so passive tracer is here, but I guess it just for show, nothing to do with twophase*/


scalar f[];
scalar * tracers = {f};
/*injection on the RIGHT (I flipped the website's description)*/

u.n[right]=dirichlet(1.);
p[right]=neumann(0.);
pf[right]=neumann(0.);
f[right]=dirichlet(y<0);

/*outflows*/
u.n[left]=neumann(0.);
p[left]=dirichlet(0.);
pf[left]=dirichlet(0.);

/*condition for the (yet nonexistant) cylinder*/
bid cylinder;
u.t[cylinder]=dirichlet(0.);

/*I think it's init*/
event init(t=0){
  mask ( y >  0.5 ? top:
        y < -0.5 ? bottom:
        sq(x)+sq(y) < sq(0.0625) ? cylinder: none);
  
#if 1 // this crashes the preprocessor...
  const face vector muc[]={0.0007825,0.00078125};
#else // ... while this does not
  const face vector muc[] = {0.0007825,0.00078125};
#endif
  
  mu=muc;
  foreach() {
    u.x[]=-1.;
    f[]=0.;
  }
  boundary({f});
}

event logfile(i++){fprintf(stderr,"%d %g %d %d \n",i,t,mgp.i,mgu.i);}

event adapt (i++) {adapt_wavelet ({u,f}, (double[]){3e-2,3e-2,3e-2}, 9, 4);}

/*domain shape*/
int main() {
  L0=8.;
  origin(-0.5,-L0/2.);
  N=512;
  run();
}
