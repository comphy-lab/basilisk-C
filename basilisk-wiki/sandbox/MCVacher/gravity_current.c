/**
# Gravity current

To be described...

An analog of the [splouch](/sandbox/M1EMN/Exemples/column_viscous.c), here with twophase.h.

Note : remove_droplets needs to be added.

## Simulation
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h" 

// I know understood that embed.h + twophase.h + "gravity" => need for reduced.h to impose gravity. 


double h_1;
double l_1;
double mu_b;
double rho_b;
double R;

// Height of the solid ceiling
double ceiling;   

// Domain extents used for visualization
double xmin, ymin;
double xmax, ymax;

FILE * fpmax; //

int main() {

  L0=4; //size of the box
  rho2 = 1000;
  rho1=rho2*1.01;
  mu1 = 0.001;
  mu2 = mu1;
  h_1=0.3; 
  l_1=0.3;
  ceiling=1.5*h_1;

  TOLERANCE = 1e-3 [*];
  
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = neumann(0.);
  
  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);
  
  u.n[left] = dirichlet(0.);
  u.t[left] = neumann(0.);
  
  u.n[right] = dirichlet(0.);
  u.t[right] = neumann(0.);
  
  G.y = -9.81;
 
  N=512;
  origin (0, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w");

  run();
}

event init (t = 0) {
  mu_b=mu2/mu1;
  rho_b=rho2/rho1;
  R=h_1/l_1;
  
  xmin = 0;
  xmax = L0;
  ymin = 0.0;
  ymax = ceiling;

  char dim_adim[80];
  sprintf (dim_adim, "dim_adim.txt"); //to compute the non-dim. parameters
  FILE * fparam = fopen(dim_adim, "w");
  fprintf (fparam, "%g %g %g\n",mu_b,rho_b,R);
  fclose (fparam);

  fraction (f, intersection((l_1-x),(h_1-y)));
  solid (cs, fs, y<ceiling);
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 0.01; t <= 30){

  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

/**
We generate a video of the fraction ...
*/

event ppm_output (t = 0; t += 0.01; t <= 30) {

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  //u.y
  char name_u[80];
  sprintf (name_u, "u_y.mp4");
  output_ppm (u.y, file = name_u, n = 512,min = -0.3, max = 0.3, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  //u.y
  char name_u_2[80];
  sprintf (name_u_2, "u_x.mp4");
  output_ppm (u.x, file = name_u_2, n = 512,min = -0.3, max = 0.3, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
}

/**
![Gravity current](gravity_current/f.mp4)
![Vertical velocity](gravity_current/u_y.mp4)
![Horizontal velocity](gravity_current/u_x.mp4)
*/

/**
## Comparison with theory : ...


*/


