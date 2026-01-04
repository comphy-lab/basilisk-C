/**
# Unstable 2D Jet.

Schlichting's 2D jet (see [here](/sandbox/MCVacher/schlichting.c)) is unstable after $Re \approx 30$ if you allow the left-right symmetry to be broken. 

It seems that it is a Kelvin-Helmotz instability created by the shear on the sides of the jet that is the physical mechanism behind it. Starting from a fluid at rest, this shear instability can be created by the initial bulb of the jet when issuing the nozzle. Here you can admire the beautiful vortices it creates...

*/

//test schlichting

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; //

face vector muv[];

int main() {

  Re=40;

  R_d=0.005; 
  L0=1; 
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);

  s[bottom] = dirichlet (U0*(x > -R_d && x <R_d));

  u.n[left] = neumann(0);
  p[left] = dirichlet(0);
  
  u.n[right] = neumann(0);
  p[right] = dirichlet(0);

  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  char param_dim[80];
  sprintf (param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf (fparam, "%g %g %g %d\n",L0,R_d,U0,N);
  fclose (fparam);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
}

event logfile (i++) { 
  fprintf (stderr, "%d %g\n", i, t); 
  fprintf (fpmax, "%d %g\n", i, t);
}

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g %g %g", x, y, u.x[], u.y[], p[],s[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}


event ppm_output (t = 0; t += 0.05; t <= 50) {

  char name2[80];
  sprintf (name2, "uY.mp4");
  output_ppm (u.y, file = name2, n = 512, min = -U0, max = +U0, linear = true);

  // optionally tracer
  char name3[80];
  sprintf (name3, "s.mp4");
  output_ppm (s, file = name3, n = 512, min = 0., max = U0, linear = true);
}

/**
![Passive tracer](unstable_2D_jet/s.mp4)
![Vertical velocity](unstable_2D_jet/uY.mp4)
*/