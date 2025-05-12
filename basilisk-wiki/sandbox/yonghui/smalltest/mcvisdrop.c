/**
##water drop in viscous environment by momentum conservation methodes
Based on [rising.c](http://mail.basilisk.fr/src/test/rising.c),

We replace two liquids with each other in order to simulate falling drops of water in a viscous environment.

The main part of the code has not been modified, comparisons will be neglected here. We simply assume the results is good.
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

u.t[right] = dirichlet(0.);
u.t[left]  = dirichlet(0.);
u.n[right] = dirichlet(1.45);
u.n[left]  = dirichlet(1.45);

p[right] = dirichlet(0);
p[left]  = neumann(0);
pf[right] = dirichlet(0);
pf[left]  = neumann(0);

uf.n[bottom] = 0.;
uf.n[top] =0.;
/**
## Main 
We need to set the calculation domain large enough to avoid the drop contact boundaries
*/
int main() {
  size (1.5);
  init_grid (128);
  rho1 = 1000., mu1 = 10.;
  rho2 = 100., mu2 = 1., f.sigma = 24.5;
  TOLERANCE = 1e-4;
  run();
}
/**
## initial & acceleration
*/
event init (t = 0) {
  mask (y > 0.6 ? top : none);
  refine (sq(x - 0.5) + sq(y) - sq(0.25) < 0. && level < LEVEL);
  fraction (f, -1.*(sq(x - 0.5) + sq(y) - sq(0.25)));
  output_ppm (f, linear = true, n=512, file = "init.png");	
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
  boundary({p});
}

/**
## some output event. 
I upload a temporary results, as for now , we get a same final shape as the one with level = 12 & solved by [visdrop.c]  

For the moment, it's good, I'll upload the necessary studies later
*/

event outputvideo (i += 20 ) {
  output_ppm (f, linear = true, n=1024,file = "f.mp4");
  output_ppm (u.x, linear = true, n=1024, file = "ux.mp4");
  output_ppm (p, linear = true, n=1024, file = "p.mp4");
  fprintf(stderr,"%d %g\n", i,t);
}

event interface (t = 7) {
  output_facets (f, stderr);
}

int k=0;
event velo (t += 0.5) {
  if (t > 2.){
    k += 1;
    char name[40];
    char velofield[40];
    sprintf (velofield, "ve-%d", k);
    sprintf (name, "dump-%d", k);
    dump (file = name);
    output_field((scalar *){u}, fopen (velofield, "w"));
  }
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL,6);
}

/**
I will update the one can be done by ploting code later

![Compare](/sandbox/yonghui/smalltest/comp1.png)
*/

