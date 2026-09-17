/**
# Self-Induced Sloshing by a jet with a cross-current

Similar to [slosh_reg.c](https://basilisk.fr/sandbox/MCVacher/slosh_reg.c) but with a global cross-current going from middle left to bottom right.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tag.h"

double h;
double U0;
double R_d;
double grav;

FILE * fpmax; //

int main() {

  R_d=0.003; 
  L0=0.4; 
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.8;
  h=0.15;
  grav=9.81;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
  p[top] = dirichlet(0.);
  
  // for left boundary, a vertical jet is added to have a global flow from left to right. Global volume of fluid is constant on average.

  u.n[left] = (y >= 0. && y <= R_d) ? dirichlet(-U0) : (y >= 2.*h/3. - R_d && y <= 2.*h/3. + R_d) ? dirichlet(0.5*U0) : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(2*U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);
 
  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w");
  
  f.sigma = 0.072;

  run();
}

event init (t = 0) {
  fraction (f, y<h);
}

/**
We initiate gravity, which is opposing to the inertia of the jet:
*/

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
We save interfaces for complex orthogonal decomposition (to find the solshing modes). We can also track at each height $y$ the maximum of $|\underline{u}|$ and have the "position" of the jet through time, and apply the same post-treatment. A version with harmonics.h should be uploaded soon...
*/

int isave1 = 1;
event res_save (t += 0.05; t <= 20) {
  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);
  
  isave1++;
}

/**
We generate videos:
*/
event ppm_output (t = 0; t += 0.05; t <= 20) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf (name2, "uX.mp4");
  output_ppm (u.x, file = name2, n = 512, min = -U0, max = U0, linear = true);
  
  scalar omega[];
  vorticity (u, omega);
  
  char name3[80];
  sprintf (name3, "omega.mp4");
  output_ppm (omega, file = name3, n = 512, linear = true);
  
}

/**
![Free-surface](slosh_reg_cross_current/f.mp4)
![Vertical velocity](slosh_reg_cross_current/uY.mp4)
![Horizontal velocity](slosh_reg_cross_current/uX.mp4)
![Vorticity](slosh_reg_cross_current/omega.mp4)
*/


