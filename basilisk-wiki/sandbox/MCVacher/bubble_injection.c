/**
# Two-Phase Bubble Plume Injection into a Liquid Bath with Embedded Geometry

This simulation extends the single-phase jet sloshing configuration of
[slosh_reg](https://basilisk.fr/sandbox/MCVacher/slosh_reg.c) to a
two-phase framework. Instead of a liquid jet impinging on a free surface, 
a gas bubble plume is injected upward from a nozzle at the bottom of a 
liquid bath.

The geometry includes two solid walls on either side of the nozzle inlet to generate bubbles continuously. If better methods exists, I am glad to be taught !
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h" 
#include "reduced.h"
#include "navier-stokes/perfs.h"

double h;
double h_2;
double U0;
double R_d;
double grav;

FILE * fpmax; 

int main() {

  R_d=0.005; // Jet nozzle radius
  L0=0.4;   
  rho2 = 10.;
  rho1=1000.;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.2;
  h=0.12;
  h_2=0.05;
  grav=9.81;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet ((x > -R_d && x <R_d)*U0*(1-(x/R_d)*(x/R_d)));
  u.t[bottom] = dirichlet(0.);
  f[bottom] = !((x > -R_d && x <R_d));

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);
  pf[top] = dirichlet(0.);

  u.n[left] = dirichlet(0.);
  u.n[right] = dirichlet(0.);

  u.n[embed]=dirichlet(0.);
  u.t[embed]=dirichlet(0.);
 
  G.y = -grav;
 
  N=128;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w"); 

  f.sigma=0.2;

  run();
}

event init (t = 0) {
  fraction (f, union(union(intersection((y<h+h_2),(y>h_2)),intersection((x<-R_d),(y<h_2))),intersection((x>R_d),(y<h_2))));
  solid(cs,fs, !union(intersection((x<-R_d),(y<h_2)),intersection((x>R_d),(y<h_2))));
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
We generate videos:
*/
event ppm_output (t = 0; t += 0.02; t <= 10) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, linear = true);
  
  char name2[80];
  sprintf (name2, "cs.mp4");
  output_ppm (cs, file = "cs.mp4", n = 512, min = 0, max = 1, linear = false);
}

/**
![Free-surface](bubble_injection/f.mp4)
![Vertical velocity](bubble_injection/uY.mp4)
![Geometry](bubble_injection/cs.mp4)
*/