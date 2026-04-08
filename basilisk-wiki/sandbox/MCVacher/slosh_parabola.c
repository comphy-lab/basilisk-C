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
#include "tracer.h"
#include "harmonic.h"

double omegas[3];

double h;
double h_1;

FILE * fpmax; //

int main() {

  L0=0.4; //size of the box
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.5;
  mu2 = 0.001*mu1;
  h=0.15; //water height at rest
  
  h_1=0.02; //size of the perturbation

  TOLERANCE = 1e-3 [*];
  
  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);

  G.y = -9.81;
  
  N=128;
  origin (-L0/2, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w");

  run();
}

event init (t = 0) {

  fraction (f, h_1*sin(3*M_PI*x/L0)+h-y);
  const face vector G[] = {0,-grav};
  
  /**
  We compute the theoretical sloshing frequencies (pseudo-pulsation in case of viscosity)
  */
  
  double lamb1 = 2.0*mu1/rho1*pow(M_PI/L0,2);
  double lamb3 = 2.0*mu1/rho1*pow(3.0*M_PI/L0,2);
  double f1 = sqrt(1.0*grav/(4.0*M_PI*L0) * tanh(M_PI*h/L0));
  double f3 = sqrt(3.0*grav/(4.0*M_PI*L0) * tanh(3.0*M_PI*h/L0));
  omegas[0] = 2.0*M_PI*f1*(sqrt(1-pow(lamb1/(2.0*M_PI*f1),2)));
  omegas[1] = 2.0*M_PI*f3*(sqrt(1-pow(lamb3/(2.0*M_PI*f3),2)));
  //omegas[0] = 2.0*M_PI*f1;
  //omegas[1] = 2.0*M_PI*f3;
  omegas[2] = 0.0;
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}


/**
 We use harmonic.h and perform it on u.y at the chosen $\omega_n$.
*/

event harmonic_analysis (t += 0.01; t <= 15) {
  harmonic_decomposition(u.y, t, omegas);
}


event profile (t = end) {
  printf ("-----END-----\n");
}

/**
We track the interface to compute the logarithmic decrease linked to the viscous dissipation.
*/

int isave1 = 1;
event res_save (t += 0.005; t <= 15){

  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

/**
We generate a video of the fraction to visualize the sloshing, as well as videos using "harmonic.h".
*/

event ppm_output (t = 0; t += 0.01; t <= 15) {

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);

  // amplitude for mode 1 (omegas[0])
  scalar amp_mode_1[];
  int k1 = 0;
  foreach()
    amp_mode_1[] = sqrt(sq(val(u.y.harmonic.A[k1])) + sq(val(u.y.harmonic.B[k1])));
  boundary ({amp_mode_1});
  char name_u_1[80];
  sprintf(name_u_1, "u_y_1.mp4");
  output_ppm(amp_mode_1, file = name_u_1, n = 512, min = 0, max = 0.1, linear = true);

  // amplitude for mode 3 (omegas[1])
  scalar amp_mode_3[];
  int k3 = 1;
  foreach()
    amp_mode_3[] = sqrt(sq(val(u.y.harmonic.A[k3])) + sq(val(u.y.harmonic.B[k3])));
  boundary ({amp_mode_3});
  char name_u_3[80];
  sprintf(name_u_3, "u_y_3.mp4");
  output_ppm(amp_mode_3, file = name_u_3, n = 512, min = 0, max = 0.1, linear = true);
  
  //u.y
  char name_u[80];
  sprintf (name_u, "u_y.mp4");
  output_ppm (u.y, file = name_u, n = 512,min = -0.3, max = 0.3, linear = true);
}

/**
![Waves](slosh_dissipation/f.mp4)
*/

/**
![Vitesse verticale](slosh_parabola/u_y.mp4)
![Projection de la norme de la vitesse verticale sur le mode n°3](slosh_parabola/u_y_3.mp4)
![Projection de la norme de la vitesse verticale sur le mode n°1](slosh_parabola/u_y_1.mp4)
*/