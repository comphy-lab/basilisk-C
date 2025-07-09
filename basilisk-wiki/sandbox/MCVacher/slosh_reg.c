/**
##Self-Induced Sloshing by a Jet

It is an adaptation from a similar experiment to 'Mechanism of jet-flutter : Self-Induced Oscillation of an Upward Plane Jet Impinging on a Free-Surface', *Madarame et al.*, 1998. Viscosity is 100 times higher than for water, showing that the Re has little influence on the instability threshold ($Fr \approx 0.6$). Grid is regular.

The code has dimensions as it is right now.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "navier-stokes/conserving.h"
#include "vtknew.h" //paraview visualisation (from Fulster's sandbox)
#include "tracer.h"

double h;
double U0;
double R_d;

/**
A passive tracer is injected with the jet to follow lagrangian paths.
*/

scalar s[];
scalar * tracers = {s};

int main() {

  R_d=0.005; 
  L0=0.4; 
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.6;
  h=0.15;

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
 
  N=64;
  origin (-L0/2, 0); //set the origin
  init_grid(N);

  run();
}

/**
We initiate gravity, which is opposing to the inertia of the jet.
*/

event init (t = 0) {
  fraction (f, y<h);
  const face vector G[] = {0,-9.81};
  a=G;
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
We save interfaces and all the vector fields in .txt files (I process data with *python* after for now, but some embedded Basilisk direct post-treatments will be added later on).
*/

int ivtk = 1;
event res_save (t += 0.01; t <= 100) 
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  
  sprintf (name, "interface-%d.txt", ivtk);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  sprintf (name, "res-%d.txt", ivtk);
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[], f[],s[],omega[]);
  fclose(fpfacet);

  ivtk++;
}

/**
.vtk snapshots are generated 10 times per seconds to visualize the velocities in real time. 
*/

int ivtk2 = 1;
event movies (t += 0.1; t <= 100) 
{
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  sprintf (name, "snapshot-%d.vtk", ivtk2);
  FILE * fpvtk = fopen(name, "w");
  output_vtk ({omega,u.x,u.y,p,f,s}, fpvtk);
  fclose(fpvtk);

  ivtk2++;
}

/**
#Video : Passive tracer

[Click me](https://drive.google.com/file/d/1RvgZekWZUlp5pdYI7DVEN7V2ZgJgbOpr/view?usp=drive_link)

<!--[![try](https://markdown-videos-api.jorgenkh.no/url?url=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3DcfL_EcfYiFM)](https://www.youtube.com/watch?v=cfL_EcfYiFM)-->

*/

