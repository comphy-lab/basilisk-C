/**
# 1D stationnary Diffusion test. */

//Geometric parameters of the drop
#define H (0.3)

// Physical parameters
#define Dr 1.
#define D1 1e-3
#define D2 (D1/Dr)

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

#define MAXTIME (30)

scalar c[];

c[left] = dirichlet(1);
c[right]  = dirichlet(0);

int maxlevel;
int minlevel = 7;
double dt;
double initial_c_mass;
mgstats mgd1;

int main ()
{
  size(H);
  origin (0.);
  DT = 0.0036; 
  TOLERANCE = 1e-4;
  init_grid (1 << minlevel);
  run();
}

event init (i = 0) {
  foreach(){
    c[]=((x<=H/2.)? 1 : 0);
  }
  boundary ((scalar *){c});
  initial_c_mass = statsf(c).sum;
}

event time_integration (i++) {
  dt = dtnext (DT);
  face vector D[];
  foreach_face(x)
    D.x[]=((x<=H/2.)? fm.x[]*D2 : fm.x[]*D1);
  boundary ((scalar *){D});
  mgd1 = diffusion (c, dt, D); 
}

event perfs (i++) {
  if (i == 0)
    fprintf (ferr, "t dt mgd1.i grid->tn perf.t perf.speed npe\n");
  fprintf (ferr, "%g %g %d %ld %g %g %d\n", 
      t, dt, mgd1.i, grid->tn, perf.t, perf.speed, npe());
}

//event adapt (i++) {
//  adapt_wavelet ({c}, (double[]){0.01}, maxlevel, minlevel);
//}

event outputs (t +=0.1; t <= MAXTIME) { 
  double total_c = statsf(c).sum;
  static  FILE * fp1 = fopen ("cmass.dat", "w");
  fprintf (fp1, "%g %g %g %g\n", t, (total_c - initial_c_mass), D1, D2);
  fflush(fp1);
}

/**

## Results
~~~gnuplot Variation of tracer quantity from the initial condition.
reset
set format y '%.2e'
D2=0.001
e=0.3
set xlabel 't* [-]'
set ylabel 'C-C0'
plot 'cmass.dat' u ($1*D2/(e**2)):2 w l t ''
~~~
 */
