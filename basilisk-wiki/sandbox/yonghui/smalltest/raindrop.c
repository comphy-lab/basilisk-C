/**
##raindrop test (Bugs still exist)
We want to simulate the shape of the raindrops. To get results faster, we use the momentum conservation method and the axisymmetric assumption.
*/
#include "axi.h" // fixme: does not run with -catch
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
# define LEVEL 8
# define UIN 3.48

/** 
## Boundary conditions */

u.n[left]  = dirichlet(UIN);
p[left]   = neumann(0.);
pf[left]  = neumann(0.);

u.n[right] = neumann(0.);
p[right]    = dirichlet(0.);
pf[right]   = dirichlet(0.);

uf.n[bottom] = 0.;
uf.n[top] = 0.;
/** 
## Main loop */
int main() {
  size (0.01);
  init_grid (64);

  rho1 = 1000.;
  rho2 = 1.2;

  mu1 = 0.001;
  mu2 = 2e-5;
  f.sigma = 0.073;

  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
  mask(y > 0.002? top :none);
  /**  The bubble is centered and has a radius of 2mm */
  refine ( sq(x - L0/2.)  + sq(y) - sq(0.0011) < 0 &&
          sq(x-L0/2.)  + sq(y) - sq(0.0009) > 0 && level < 8);
  fraction (f, -1.*(sq(x - L0/2.) + sq(y) - sq(0.001)));
  foreach()
    u.x[] = UIN*(1.-f[]);
  output_ppm (u.x, linear = true,min = -2*UIN, max = 2*UIN, n=512, file = "uxinit.png");	
}



/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 9.8;
}

/**
We output the shape of the bubble every iteration. */
event animation ( i +=10; t <= 0.1){
  output_ppm (f, linear = true, n=512, file = "f.mp4");
  output_ppm (u.x, linear = true,min = -10, max = 10, n=512, file = "ux.mp4");
  output_ppm (p, linear = true, n=512, file = "p.mp4");
  fprintf (stderr, "%d %g %g \n", i, t, dt);
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-1,1e-1}, LEVEL,4);
}
