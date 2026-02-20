#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "tracer.h"
//#include "tag.h"

double h;
double U0;
double R_d;

double t_period;
double t_max;

double Fr;
double Re;
double mu_b;
double rho_b;
double R_1;
double R_2;

scalar s[];
scalar * tracers = {s};

FILE * fpmax; 

int main() {

  R_d=0.004; //Rayon du jet initial
  L0=0.8; //Taille de la boite  
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.75;
  h=0.1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  p[bottom] = dirichlet(0.);
  pf[bottom] = dirichlet(0.);

  u.n[left] = neumann(0.);
  u.n[top] = neumann(0.);
  u.n[right] = neumann(0.);

  N=256;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w"); 

  f.sigma = 0.072;
  G.y=-9.81;

  run();
}

event init (t = 0) {
  fraction (f, y<h);
  boundary ({f});
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event ppm_output (t = 0; t += 0.01; t <= 30) { //t_max
  
char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);

  char name4[80];
  sprintf (name4, "uX.mp4");
  output_ppm (u.x, file = name4, n = 512, min = -U0, max = +U0, linear = true);
  
  char name5[80];
  sprintf (name5, "p.mp4");
  output_ppm (p, file = name5, n = 512, min = 0, max = 0.1*U0*U0, linear = true);
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
![Pressure](jet_flutter/p.mp4)
*/