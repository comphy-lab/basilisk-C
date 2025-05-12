/**
# Rayleigh-Taylor instability (homework given by F.Hecht)
Used to comapre with Freefem++ code, it's a very simple version since im a newbee in Freefem++
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#define tend 3.

int main() {
  N = 256;
  periodic(right);
  rho1 = 1., mu1 = 0.0001;
  rho2 = 2., mu2 = 0.0001;
  f.sigma=0.00001; // put any value you want
  run(); 
}

/**
## boundary conditions */

u.n[left]  = dirichlet(0.);
u.n[right] = dirichlet(0.);

u.n[top]  = dirichlet(0.);
u.t[top]  = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);


p[bottom]    = neumann(0.);
pf[bottom]   = neumann(0.);

p[top]   = dirichlet(0.);
pf[top]  = dirichlet(0.);

/**
## init event */

event init (t = 0) {
  DT = 0.01;
  fraction (f,  0.5 + 0.02*cos(x*2.*M_PI) - y );
  output_ppm (f, file = "finit.png", linear = true, n=512, min = 0, max = 1);
}


event movies (i += 1; t <= tend) {
  /**  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "vort.gif",min = -10, max = 10, linear = true);
*/
  output_ppm (f, file = "f.mp4", linear = true, min = 0, max = 1);
  fprintf(stderr,"%d %g %g\n" , i, t, dt);
}

event acceleration (i++; t <= tend) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 0.98;
  //a possible change to RM instability
  //av.y[] += ( y <= (t + 0.05) && y >= ( t)  ) ? 25.6 :0. ;  
}

 
We adapt according to the error on the velocity and tracer fields. 
event adapt (i++) {
  adapt_wavelet ({u,f}, (double[]){3e-2,3e-2,3e-2}, 8, 4);
} 