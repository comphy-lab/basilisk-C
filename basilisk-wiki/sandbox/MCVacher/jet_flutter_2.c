#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tracer.h"
#include "reduced.h"
#include "tag.h"

double h;
double U0;
double R_d;

scalar s[];
scalar * tracers = {s};

FILE * fpmax; 

#define MAXLEVEL 8

int main() {

  R_d=0.02; //Rayon du jet initial
  L0=2; //Taille de la boite  
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.12;
  mu2 = 0.001*mu1;
  U0=0.8;
  h=0.3;

  TOLERANCE = 1e-3 [*];
  
  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.n[top]= neumann(0.);

  u.n[left] =  y < h ? dirichlet(0.) : neumann(0.); 
  u.n[right] = y < h ? dirichlet(0.) : neumann(0.);
  
  N=128;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w"); 

  f.sigma = 0.072;
  G.y= -9.81;
  Z.y = h;

  run();
}

event init (t = 0) {
  fraction(f, h-y);
}

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event ppm_output (t = 0; t += 0.01; t <= 15) {
  output_ppm (f, file = "f.mp4", n = 512, min = 0, max = 1, linear = true);
  output_ppm (u.y, file = "uY.mp4", n = 512, min = -U0, max = +U0, linear = true);
  output_ppm (s, file = "s.mp4", n = 512, min = 0., max = U0, linear = true);
  output_ppm (u.x, file = "uX.mp4", n = 512, min = -U0, max = +U0, linear = true);
  output_ppm (p, file = "p.mp4", n = 512, linear = true);
}

event adapt (i++){
  adapt_wavelet((scalar*){f,u}, (double[]){0.01,0.05,0.05}, MAXLEVEL,5);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
# Outputs
![Free-surface](jet_flutter_2/f.mp4)
![Passive tracer](jet_flutter_2/s.mp4)
![Vertical velocity](jet_flutter_2/uY.mp4)
![Horizontal velocity](jet_flutter_2/uX.mp4)
![Pressure](jet_flutter_2/p.mp4)
*/