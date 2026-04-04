/**
# Low Reynolds jet confined in a box (bottom aspiration)

This simulation models a 2D planar jet confined inside a rectangular box. A low Reynolds instability of the jet takes place, reaching a limit cycle.
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
  Re      = 26;      // Reynolds number
  
  R_d     = 0.012;   // Jet radius
  L0      = 1.;      // Box size
  U0      = 1.;      // Inlet velocity
  ceiling = 0.8;     // Ceiling height

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

  xmin = -L0/2;
  xmax = L0/2;
  ymin = 0.0;
  ymax = ceiling;

  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
  }

  solid (cs, fs, y<ceiling); // Embedded ceiling using embed.h
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
Movie outputs for velocity, tracer, and pressure fields.
*/
event ppm_output (t = 0; t += 0.05; t <= 100) {
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
  output_ppm (p, file = name3, n = 512, min = 0, max = 0.1*U0*U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
}


/**
Periodic dump for restart and post-processing.
*/
double t_prev_dump = 0.;

event dump_state (t = 0; t += 10; t <= 100) { //t_max
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
  printf ("----END----\n");
}

/**
## Movies

![Passive tracer](box_jet/s.mp4)
![Vertical velocity](box_jet/uY.mp4)
![Horizontal velocity](box_jet/uX.mp4)
![Pressure](box_jet/p.mp4)
*/