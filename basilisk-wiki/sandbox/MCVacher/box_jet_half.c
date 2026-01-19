//Jet impinging on a wall (see Bouchet and Climent 2012)

//L=10D = 20R, exit=2d=4R

//for now, it is just a half-domain (trying to get base flow without the instability

//only pressure of the inside of the box is driving the flow outside on the left

//right now the boundary condition are not really good

//Notes : implement parabolic inlet, compare with the paper etc...

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; //

face vector muv[];

int main() {

  Re=250;

  R_d=0.5; 
  L0=10; 
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (U0*(x > -R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = dirichlet(0);
  u.t[top] = dirichlet(0);

  s[bottom] = dirichlet (U0*(x > -R_d));
  p[bottom] = dirichlet(0.);


  u.n[left] = (y<4*R_d && u.n[] < 0.) ? neumann(0) : dirichlet(0);

  u.n[right] = dirichlet(0.);
  p[right]=neumann(0.);

  N=128;  
  origin (-L0, 0);
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
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t); 
  fprintf (fpmax, "%d %g \n", i, t);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

event ppm_output (t = 0; t += 0.5; t <= 400) {

  char name2[80];
  sprintf (name2, "uY.mp4");
  output_ppm (u.y, file = name2, n = 512, min = -U0, max = +U0, linear = true);

  char name3[80];
  sprintf (name3, "s.mp4");
  output_ppm (s, file = name3, n = 512, min = 0., max = U0, linear = true);
}

/**
![Passive tracer](box_jet_half/s.mp4)
![Vertical velocity](box_jet_half/uY.mp4)
*/