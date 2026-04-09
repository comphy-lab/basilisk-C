/**
# Sloshing in a parabolic container (NS)

DNS of the Saint-Venant version of [parabola.c](https://basilisk.fr/src/test/parabola.c)

## Simulation
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tag.h"
//#include "harmonic.h"

double h;
double h_1;
double amp;

FILE * fpmax; //

int main() {

  L0=1; //size of the box
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.05;
  mu2 = 0*mu1;
  h=0.31; //water height at rest
  amp=0.8;
  
  h_1=0.12; //size of the perturbation

  TOLERANCE = 1e-3 [*];
  
  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);

  G.y = -9.81;
  
  N=256;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w");

  run();
}

event init (t = 0) {
  fraction (f, intersection(h_1*x+h-y,-(x/amp)*(x/amp)+y-0.02));
  solid (cs, fs, -(x/amp)*(x/amp)+y-0.02);
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 0.02; t <= 10){

  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

event ppm_output (t = 0; t += 0.02; t <= 10) {

  //Note for me : do some tracers only for the water
  
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name_u_x[80];
  sprintf (name_u_x, "u_x.mp4");
  output_ppm (u.x, file = name_u_x, n = 512, min = -0.3, max = 0.3, linear = true);
  
  char name_u_y[80];
  sprintf (name_u_y, "u_y.mp4");
  output_ppm (u.y, file = name_u_y, n = 512,min = -0.3, max = 0.3, linear = true);
}

/**
![Waves](slosh_parabola/f.mp4)
![Vertical velocity](slosh_parabola/u_y.mp4)
![Horizontal velocity](slosh_parabola/u_x.mp4)
*/