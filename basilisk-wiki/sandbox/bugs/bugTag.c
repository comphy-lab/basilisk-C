/**

This Code produce a strange behaviour in the tag function. The function didn't
mark correctly the fluid. To observe correctly the bug, you need to use gfsview.

The interresting parameter to observe is m. m should only change from 0 to 
something and from something to 0.*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

/**
The bug was observe for a high level of refinement in the original case.*/

#define LEVEL 13
#define L0 10
#define R1 1.
#define L2 0.05
#define L1 (sqrt(sq(R1)-4*sq(L2)) - R1 - 2*L2)
#define Sigma 1.

double uemax = 0.35;

int main(int argc, char *argv[]) {
  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << (5));


  double La = 500.;
  
  double Oh = sqrt(1./La);


  rho1 = 1., mu1 = Oh;
  rho2 = 1./998., mu2 = Oh/55;

  f.sigma = Sigma;
  
  run();
}

/**
This geometry is linked to the bug. We observe the bug with the geometry of 
the code [bursting.c](/sandbox/aberny/bursting.c). We just reproduce this 
geometry here.*/

double geometry(double x, double y) {
  
  double C1 = sq(x + R1 + 2*L2) + sq(y) - sq(R1);
  double C2 = sq(x + L2+0.003125) + sq(y - 2*L2) - sq(L2+0.003125);
  double D1 = - x - 1e-8;
  double D2 = y - (2)*L2;
  double D3 = - x + (sqrt(sq(R1)-4*sq(L2))-R1-2*L2);

  double D1D2 = min(D1, D2);
  double D1D2D3 = max(D1D2, D3);
  double D1D2D3C1 = min(D1D2D3, C1);
  double D1D2D3C1C2 = min(-D1D2D3C1, C2);
  return -D1D2D3C1C2;
}


event init (t = 0) {
  double iteration = 0;
  do {
    iteration++;
    fraction(f, geometry(x,y));
  } while( adapt_wavelet({f,u}, (double []){0.01,uemax,uemax,uemax},
			 maxlevel = LEVEL, 5).nf != 0 && iteration <= 10);
  static FILE * fp = fopen("initial.png", "w");
  output_ppm(f, fp, 400, min = 0, max = 1);
}

event adapt (i ++) {
  adapt_wavelet ({f,u}, (double []){0.01,uemax,uemax,uemax},
		 maxlevel = LEVEL, 5);
}

static double compt = 0;

/**
The time 0.55 correspond to the time for the simulation to produce a jet, with a droplet. If there is no more bug, then the simulation will stop at t=0.55*/

event bug (i ++; t<=0.55) {
  scalar m[];
  foreach()
    m[] = f[] > 1.5e-1;
  int n = tag(m);
  char name[80];
  sprintf(name,"bug-%05ld.gfs",i);

  /**
  If the tag function return a value higher than 1, then we are suppose to 
  have, in the simulation, a droplet.*/

  if (n>1) {
    compt++;
    FILE* fp = fopen (name, "w");
    output_gfs(fp);
    fprintf(stderr, "droplet? %d %d %g \n", compt, i, t);
  }

  printf("%d %g %d\n", i, t, n);
  fflush(stdout);

  /**
  We just output in a gfs file 5 steps before stopping the simulation.*/
  if (compt>5) {
    return 1;
  }
}

event logfile (i++) {
  fprintf (stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgpf.i, mgu.i);
}

event finalState (t = end) {
  scalar m[];
  foreach()
    m[] = f[] > 1.5e-1;
  int n = tag(m);
  output_gfs (file = "final.gfs");
  static FILE * fp = fopen("final.ppm", "w");
  output_ppm(m, fp, 512, min=0, max=1);
  dump (file = "dump");
}
