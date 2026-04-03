/**
# Title

Bla bla
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"

// Passive tracer
scalar s[];
scalar * tracers = {s};

// Physical and geometrical parameters

double U0;        // Jet inlet velocity
double R_d;       // Jet radius
double Re;        // Reynolds number

FILE * fpmax; 

face vector muv[];

int main() {

    // Flow parameters
  Re      = 60;      // Reynolds number
  R_d     = 0.02;   // Jet radius
  L0      = 1.;      // Box size
  U0      = 1.;      // Inlet velocity
  
  TOLERANCE = 1e-3 [*];

  
  // Boundary conditions
  
  // Bottom inlet: top-hat velocity profile
  u.n[bottom] = dirichlet(U0*(x > -R_d && x < R_d));
  u.t[bottom] = dirichlet(0.);

  // Top outlet
  u.n[top] = (x > -5*R_d && x < 5*R_d) ? neumann(0.) : dirichlet(0.);
  
  // Tracer injected with the jet
  s[bottom] = dirichlet(U0*(x > -R_d && x < R_d));
  s[top] = (x > -5*R_d && x < 5*R_d) ? neumann(0.) : dirichlet(0.);

  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  
  // Grid initialization
  N=128;  
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

/**
Viscosity definition: $\mu = U_0 \times R_d / Re$ inside the fluid.
*/
event properties (i++)
{
  foreach_face() {
    muv.x[] = fm.x[]*U0*R_d/Re;
  }
}

/**
Initialization: domain extents and embedded ceiling.
*/
event init (t = 0) {

  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
  }
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
Movie outputs for velocity, tracer, and pressure fields.
*/
event ppm_output (t = 0; t += 0.1; t <= 100) {
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -U0, max = U0, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);

  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, min = 0, max = 0.1*U0*U0, linear = true);
}

event profile (t = end) {
  printf ("----END----\n");
}

/**
## Movies

![Passive tracer](box_jet_2/s.mp4)
![Vertical velocity](box_jet_2/uY.mp4)
![Horizontal velocity](box_jet_2/uX.mp4)
![Pressure](box_jet_2/p.mp4)
*/