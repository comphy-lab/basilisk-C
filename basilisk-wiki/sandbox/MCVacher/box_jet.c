#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;
double plafond;

double xmin;
double ymin;
double xmax;
double ymax;

double Re;

FILE * fpmax; //

face vector muv[];

int main() {

  Re=30;

  R_d=0.02; 
  L0=1; 
  U0=1;
  plafond=0.8;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet(U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);
  p[embed] = dirichlet(0.);
  pf[embed]=dirichlet(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  s[bottom] = dirichlet(U0*(x > -R_d && x <R_d));

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);

  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);

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
  foreach_face() {
    muv.x[] = fm.x[]*U0*R_d/Re;
  }
}

event init (t = 0) {

  xmin = -L0/2;
  xmax = L0/2;
  ymin = 0.0;
  ymax = plafond;

  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
  }

  solid (cs, fs, y<plafond);
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event ppm_output (t = 0; t += 0.05; t <= 100) {
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = 0, max = 1, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  // optionally tracer
  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true , box = {{xmin, ymin}, {xmax, ymax}});

  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, min = 0, max = 0.1*U0*U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  char name4[80];
  sprintf (name4, "uX.mp4");
  output_ppm (u.x, file = name4, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
}

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
![Passive tracer](box_jet/s.mp4)
![Vertical velocity](box_jet/uY.mp4)
![Horizontal velocity](box_jet/uX.mp4)
![Pressure](box_jet/p.mp4)
*/