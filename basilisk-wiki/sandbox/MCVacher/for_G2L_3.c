/**
# High aspect ratio (monophasic)
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"

// Passive tracer
scalar s[];
scalar * tracers = {s};

// Physical and geometrical parameters

double U0;        // Jet inlet velocity
double R_d;       // Jet radius
double ceiling;   // Height of the solid ceiling

// Domain extents used for visualization
double xmin, ymin;
double xmax, ymax;

double Re;        // Reynolds number

FILE * fpmax; 

face vector muv[];

int main() {

    // Flow parameters
  Re      = 40;      // Reynolds number
  
  R_d     = 0.003;   // Jet radius
  L0      = 0.4;      // Box size
  U0      = 0.9;      // Inlet velocity
  ceiling = 0.8*L0;     // Ceiling height

  TOLERANCE = 1e-3 [*];

  
  // Boundary conditions
  
  // Bottom inlet: top-hat velocity profile
  u.n[bottom] = dirichlet(U0*(x > -R_d && x < R_d));
  u.t[bottom] = dirichlet(0.);
  p[embed]    = dirichlet(0.);
  pf[embed]   = dirichlet(0.);

  // No-slip condition on the embedded solid ceiling
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  // Tracer injected with the jet
  s[bottom] = dirichlet(U0*(x > -R_d && x < R_d));

  // Lateral walls: no-slip walls except for the lower part where there is imposed outflow
  u.n[left]  = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.t[left]  = y < R_d ? neumann(0.)  : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.)  : dirichlet(0.);
  
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
Viscosity definition: $\nu = U_0 \times 2*R_d / Re$ inside the fluid ("mu" in Basilisk).
*/
event properties (i++)
{
  foreach_face() {
    muv.x[] = fm.x[]*U0*2*R_d/Re;
  }
}

/**
Initialization: domain extents and embedded ceiling.
*/
event init (t = 0) {

  xmin = -L0/2;
  xmax = L0/2;
  ymin = 0.0;
  ymax = ceiling;

  solid (cs, fs, y<ceiling); // Embedded ceiling using embed.h
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
Movie outputs for velocity, tracer, and pressure fields.
*/
event ppm_output (t = 0; t += 0.01; t <= 20) {
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -U0, max = U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true , box = {{xmin, ymin}, {xmax, ymax}});

  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
}

event profile (t = end) {
  printf ("----END----\n");
}

/**
## Movies

![Passive tracer](for_G2L_3/s.mp4)
![Vertical velocity](for_G2L_3/uY.mp4)
![Horizontal velocity](for_G2L_3/uX.mp4)
![Pressure](for_G2L_3/p.mp4)
*/