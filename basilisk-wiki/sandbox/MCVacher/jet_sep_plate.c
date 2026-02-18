/**
# Jet impinging on a separating plate

Description : *soon*
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"

double h;
double h_2;
double U0;
double R_d;
double R_d_2;

/**
A passive tracer is injected with the jet to follow lagrangian paths.
*/

scalar s[];
scalar * tracers = {s};

FILE * fpmax; //

int main() {

  R_d=0.005; 
  R_d_2=0.003;
  L0=0.4; 
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.6;
  h=0.15;
  h_2=0.1;

  TOLERANCE = 1e-3 [*];

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  u.n[bottom] = dirichlet (U0*(x > -R_d && x <R_d));
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);

  s[bottom] = dirichlet (U0*(x > -R_d && x <R_d));

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);
 
  G.y = -9.81;

  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w");

  run();
}

event init (t = 0) {
  fraction (f, union(union(intersection(y<h,x<-R_d_2),intersection(y<h,x>R_d_2)),y<h_2));
  solid (cs, fs, !intersection(intersection(x>-R_d_2,x<R_d_2),y>h_2));
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t); 
  fprintf (fpmax, "%d %g \n", i, t);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 1; t <= 20) {
  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);
  
  isave1++;
}

/**
To visualize both the surface and the jet, we can follow lagrangian trajectories using a passive tracer. We generate videos:
*/
event ppm_output (t = 0; t += 0.05; t <= 20) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name0[80];
  sprintf (name0, "uX.mp4");
  output_ppm (u.x, file = name0, n = 512, min = -U0, max = U0, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);

  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, min = 0, max = 100, linear = true);
}

/**
## Movies

![Free surface motion](jet_sep_plate/f.mp4)
![Passive tracer](jet_sep_plate/s.mp4)
![Vertical velocity](jet_sep_plate/uY.mp4)
![Horizontal velocity](jet_sep_plate/uX.mp4)
![Pressure](jet_sep_plate/p.mp4)
*/
