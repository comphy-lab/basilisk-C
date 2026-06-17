/**
# High-aspect ratio (diphasic)

The code has dimensions as it is right now.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tracer.h"
#include "tag.h"

double h;
double U0;
double R_d;
double grav;
double Re;

/**
A passive tracer is injected with the jet to follow lagrangian paths.
*/

scalar s[];
scalar * tracers = {s};

FILE * fpmax; //

int main() {

  Re=40;
  R_d=0.003; 
  L0=0.4; 
  U0=0.9;
  rho2 = 1.3;
  rho1=1000;
  mu1 = 2*rho1*R_d*U0/Re;
  mu2 = 0.01*mu1;
  h=0.8*L0;
  grav=9.81;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
  p[top] = dirichlet(0.);
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);
 
  N=128;
  origin (-L0/2, 0);
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w");
  
  f.sigma = 0.072;

  run();
}

/**
We initiate gravity, which is opposing to the inertia of the jet.
*/

event init (t = 0) {
  fraction (f, y<h);
}

#if 1
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -9.81;
}
#endif

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
We kill eventual numerical bubbles (not real bubbles because there is no surface tension here.
*/

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
To visualize both the surface and the jet, we can follow lagrangian trajectories using a passive tracer. We generate videos:
*/
event ppm_output (t = 0; t += 0.01; t <= 20) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  // optionally tracer
  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

/**
![Free-surface](for_G2L_2/f.mp4)
![Passive tracer](for_G2L_2/s.mp4)
![Vertical velocity](for_G2L_2/uY.mp4)
*/

