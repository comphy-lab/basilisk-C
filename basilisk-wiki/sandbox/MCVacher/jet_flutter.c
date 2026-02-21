#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "tracer.h"

double h;
double U0;
double R_d;

scalar s[];
scalar * tracers = {s};

FILE * fpmax; 

int main() {

  R_d=0.004; //Rayon du jet initial
  L0=0.6; //Taille de la boite  
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.8;
  h=0.11;

  TOLERANCE = 1e-3 [*];

  f[left]=  y < h ? dirichlet(1.) : dirichlet(0.); 
  f[right]= y < h ? dirichlet(1.) : dirichlet(0.); 
  
  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  p[bottom] = dirichlet(0.);
  pf[bottom] = dirichlet(0.);

  u.n[left] = neumann(0.);
  u.n[right] = neumann(0.);

  N=128;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w"); 

  f.sigma = 0.072;
  G.y= -9.81;

  run();
}

event init (t = 0) {
  fraction(f, h-y);
  foreach() {
    u.x[] = 0.;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary ({f,u});
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event ppm_output (t = 0; t += 0.05; t <= 50) {
  output_ppm (f, file = "f.mp4", n = 512, min = 0, max = 1, linear = true);
  output_ppm (u.y, file = "uY.mp4", n = 512, min = -U0, max = +U0, linear = true);
  output_ppm (s, file = "s.mp4", n = 512, min = 0., max = U0, linear = true);
  output_ppm (u.x, file = "uX.mp4", n = 512, min = -U0, max = +U0, linear = true);
  output_ppm (p, file = "p.mp4", n = 512, linear = true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}
/**
# Outputs
![Free-surface](jet_flutter/f.mp4)
![Passive tracer](jet_flutter/s.mp4)
![Vertical velocity](jet_flutter/uY.mp4)
![Horizontal velocity](jet_flutter/uX.mp4)
![Face Vertical velocity](jet_flutter/ufY.mp4)
![Face Horizontal velocity](jet_flutter/ufX.mp4)
![Pressure](jet_flutter/p.mp4)
*/