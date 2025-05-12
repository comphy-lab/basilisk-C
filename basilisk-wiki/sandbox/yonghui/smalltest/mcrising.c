/**
##Momentum conservation test
Edit version of [rising.c](http://mail.basilisk.fr/src/test/rising.c), 
We set the liquid1 with a final steady velocity (-0.361) in the opposite direction to calculate the stable shape of the bubble.
*/
#include "axi.h"  
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"//actually in this simple case it's useless ??
#include "tension.h"
# define LEVEL 8
/**
## BC
*/
u.n[left] = dirichlet(-0.361);
p[left] =neumann(0.);
pf[left] = neumann(0.);

u.n[right]  = dirichlet(-0.361);
p[right]  =  dirichlet(0.);
pf[right]  = dirichlet(0.);

uf.n[bottom] = 0.;
uf.n[top] = 0.;

/**
## Main 
We need to set the calculation domain large enough to avoid the bubble contact boundaries
*/
FILE * fp;
int main() {
  fp = fopen ("forme", "w");
  size (1.4);
  init_grid (64);
  rho1 = 1000., mu1 = 10.;
  rho2 = 100., mu2 = 1., f.sigma = 24.5;
  TOLERANCE = 1e-4;
  run();
}
/**
## Initial and Gravity
*/
event init (t = 0) {
  mask (y > 0.5 ? top : none);
  fraction (f, sq(x - 0.7) + sq(y) - sq(0.25));
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}

/**
## output
The final shape will converge to the shape found in [rising.c] (update later).

Since it's just a beta version, I just put the code here.
If you are interested, please compare it by yourself.
*/
event videos (i += 2){
  //	  fprintf(stderr, "%d %g %g\n", i, t, dt);
  output_ppm (f, min = 0, max = 1,n = 1024,file = "forme.mp4");
  output_ppm (u.x, min = -0.5, max = 0.5,n = 1024,file = "ux.mp4");
  output_ppm (p, n = 1024,file = "p.mp4");
}

event finalform1 (t+=0.5 ;  t <= 5.){
  output_ppm (f, min = 0, max = 1,n = 1024,file = "forme.png");
  output_facets (f, fp);
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-3,1e-3}, LEVEL,4);
}
