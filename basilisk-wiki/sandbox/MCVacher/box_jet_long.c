/**
# High Reynolds jet in an open channel impinging on a solid surface

This simulation models a 2D planar jet entering an open channel
and impinging on a solid ceiling using the embed.h. 
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
double as_r;      // Aspect ratio of the jet

// Domain extents used for visualization
double xmin, ymin;
double xmax, ymax;

double Re;    // Reynolds number

FILE * fpmax; 

face vector muv[];

int main() {

  Re=340;
  
  as_r    = 8;                 // H/(2*R_d)
  L0      = 0.8;               // Domain size
  U0      = 1.;                // Inlet velocity magnitude
  ceiling = 0.1;               // Ceiling height
  R_d     = ceiling/(2*as_r);  // Initial jet radius

  TOLERANCE = 1e-3 [*];

  // Jet inflow at the bottom boundary (top-hat profile)
  u.n[bottom] = dirichlet(U0*(x > -R_d && x < R_d));
  u.t[bottom] = dirichlet(0.);
  p[bottom]   = dirichlet(0.);
  pf[bottom]  = dirichlet(0.);

  // No-slip condition on the embedded solid
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  // Tracer injected with the jet
  s[bottom] = dirichlet(U0*(x > -R_d && x < R_d));

  // Open boundaries on the sides
  u.n[left]  = neumann(0.);
  p[left]    = neumann(0.);
  u.n[right] = neumann(0.);
  p[right]   = neumann(0.);

  // Grid initialization
  N = 256;
  origin (-L0/2, 0);
  init_grid(N);

  // Assign spatially varying viscosity
  mu = muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}

/**
Viscosity definition : $\mu = U_0 \times 2R_d / Re$ inside the fluid.
*/

event properties (i++)
{
  foreach_face() {
    muv.x[] = fm.x[]*U0*2*R_d/Re;
  }
}

event init (t = 0) {

  xmin = -L0/2;
  xmax = L0/2;
  ymin = 0.0;
  ymax = ceiling;

  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
  }

  solid (cs, fs, y<ceiling); //Use of embed
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
Movie outputs for velocity, tracer and pressure fields
*/
event ppm_output (t = 0; t += 0.01; t <= 10) { //t_max
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true , box = {{xmin, ymin}, {xmax, ymax}});
  
  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, min = -1/2*U0*U0, max = 1/2*U0*U0 , linear = true, box = {{xmin, ymin}, {xmax, ymax}});
}

/**
Periodic dump for restart and post-processing
*/
double t_prev_dump = 0.;

event dump_state (t = 0; t += 5; t <= 10) { //t_max
  dump("restart");
  t_prev_dump = t;
  fprintf(stderr, "Dumped state at t = %g\n", t);

  char name[80];
  sprintf(name, "res_t%.1f.txt", t);
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf(fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
## Movies

![Passive tracer](box_jet_long/s.mp4)
![Vertical velocity](box_jet_long/uY.mp4)
![Horizontal velocity](box_jet_long/uX.mp4)
![Pressure](box_jet_long/p.mp4)
*/