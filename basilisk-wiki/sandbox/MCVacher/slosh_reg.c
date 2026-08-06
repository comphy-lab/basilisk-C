/**
# Self-Induced Sloshing by a Jet

It is an adaptation from a similar experiment to 'Mechanism of jet-flutter : Self-Induced Oscillation of an Upward Plane Jet Impinging on a Free-Surface', *Madarame et al.*, 1998. Viscosity is 100 times higher than for water, showing that the Re has little influence on the instability threshold ($Fr \approx 0.6$). Grid is regular.

The code has dimensions as it is right now.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tracer.h"
#include "tag.h"
#include "harmonic.h"

double h;
double U0;
double R_d;
double grav;

double omegas[3]; //to use harmonic.h

/**
A passive tracer is injected with the jet to follow lagrangian paths.
*/

scalar s[];
scalar * tracers = {s};

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
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);
 
  N=256;
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
  
  double fs1 = sqrt(1.0*grav/(4.0*M_PI*L0) * tanh(M_PI*h/L0));
  double fs3 = sqrt(3.0*grav/(4.0*M_PI*L0) * tanh(3.0*M_PI*h/L0));
  double fj1 = 1/(2*M_PI)*sqrt(grav/h);
  omegas[0] = 2.0*M_PI*fs1;
  omegas[1] = 2.0*M_PI*fs3;
  omegas[2] = 2.0*M_PI*fj1;
  omegas[3] = 0.0;
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

event harmonic_analysis (t += 0.05; t <= 70) {
  harmonic_decomposition(u.y, t, omegas);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
We save interfaces for complex orthogonal decomposition (to find the solshing modes). We can also track at each height $y$ the maximum of $|\underline{u}|$ and have the "position" of the jet through time, and apply the same post-treatment. A version with harmonics.h should be uploaded soon...
*/

int isave1 = 1;
event res_save (t += 0.05; t <= 70) {
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
event ppm_output (t = 0; t += 0.05; t <= 70) {
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
  
   // amplitude for mode 1 sloshing (omegas[0])
  scalar amp_mode_1[];
  int k1 = 0;
  foreach()
    amp_mode_1[] = sqrt(sq(val(u.y.harmonic.A[k1])) + sq(val(u.y.harmonic.B[k1])));
  boundary ({amp_mode_1});
  char name_u_1[80];
  sprintf(name_u_1, "u_y_1.mp4");
  output_ppm(amp_mode_1, file = name_u_1, n = 512, min = -0.1, max = 0.1, linear = true);

  // amplitude for jet mode (omegas[1])
  scalar amp_mode_jet[];
  int k3 = 2;
  foreach()
    amp_mode_jet[] = sqrt(sq(val(u.y.harmonic.A[k3])) + sq(val(u.y.harmonic.B[k3])));
  boundary ({amp_mode_jet});
  char name_u_3[80];
  sprintf(name_u_3, "u_y_jet.mp4");
  output_ppm(amp_mode_jet, file = name_u_3, n = 512, min = -0.1, max = 0.1, linear = true);
}

/**
![Free-surface](slosh_reg/f.mp4)
![Passive tracer](slosh_reg/s.mp4)
![Vertical velocity](slosh_reg/uY.mp4)
![Projection of vertical velocity on sloshing freq. - mode n°1](slosh_reg/u_y_1.mp4)
![Projection of vertical velocity on jet freq.](slosh_reg/u_y_jet.mp4)
*/



