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
double grav;

FILE * fpmax; //

face vector muv[];

int main() {

  Re=80;

  R_d=0.005; //Rayon du jet initial
  L0=0.5; //Taille de la boite  
  U0=0.7;
  grav=9.81;
  plafond=0.1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet(U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);
  p[bottom]=dirichlet(0.);

  u.n[embed] = dirichlet(0);
  u.t[embed] = neumann(0);

  s[bottom] = dirichlet(U0*(x > -R_d && x <R_d));

  u.n[left] = neumann(0.);
  p[left] = neumann(0.);

  u.n[right] = neumann(0.);
  p[right] = neumann(0.);


  N=256;  
  origin (-L0/2, 0);
  init_grid(N);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
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

  const face vector G[] = {0,-grav};
  a=G;
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event ppm_output (t = 0; t += 0.02; t <= 60) { //t_max
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  // optionally tracer
  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true , box = {{xmin, ymin}, {xmax, ymax}});
}

double t_prev_dump = 0.;

event dump_state (t = 0; t += 10; t <= 60) { //t_max
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
![Passive tracer](box_jet_long/s.mp4)
![Vertical velocity](box_jet_long/uY.mp4)
![Horizontal velocity](box_jet_long/uX.mp4)
*/