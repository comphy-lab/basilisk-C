/**
# 1D unstationnary Diffusion test. */

//Geometric parameters of the drop
#define H (0.3)

// Physical parameters
#define Sctw (1)
#define NUw (1e-6)
#define Dtw (1e-6)

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

#define MAXTIME (60)

scalar c[], cold[];

c[left] = neumann(0);
c[right]  = neumann(0);
  
int maxlevel;
int minlevel = 6;
double dt,errmax;
mgstats mgd1;

double cexact(double x, double t){
return (erfc((x-H/2.)/(2.*pow(Dtw*t,0.5)))/2.);
}

int main ()
{
  size(H);
  origin (0.);
  DT = 0.05; 
  TOLERANCE = 1e-4;
  for (maxlevel = 8; maxlevel <= 10; maxlevel++) {
    init_grid (1 << minlevel);
    run();
    fprintf (fout, "\n \n");
    fflush (fout);
    fprintf (ferr, "\n \n");
    fflush (ferr);
  }
}

event init (i = 0) {
  foreach(){
    c[]=((x<=H/2.)? 1 : 0);
    cold[]= c[];
  }
  boundary ((scalar *){c,cold});
}

event time_integration (i++) {
 dt = dtnext (DT);
 const face vector D[] = {Dtw};
  mgd1 = diffusion (c, dt, D); 
}

event perfs (i++) {
  if (i == 0)
    fprintf (ferr, "t dt grid->tn perf.t perf.speed npe\n");
  fprintf (ferr, "%g %g %d %ld %g %g %d\n", 
      t, dt, mgd1.i, grid->tn, perf.t, perf.speed, npe());
}

event adapt (i++) {
  adapt_wavelet ({c}, (double[]){0.01}, maxlevel, minlevel);
}

event end (t = {1,MAXTIME}) {
  scalar dce[];
  foreach()
    dce[] = c[]-cold[];
  boundary({dce});
  output_ppm (c, file = "c.png", spread=2, linear = true, n=512);
  output_ppm (dce, file = "diference_c.png", spread=2, linear = true, n=512);
  foreach()
    fprintf (fout, "%g %d %g %g\n", x, N, c[], cexact(x,t));
  fprintf (fout, "\n \n");
  fflush (fout);
}

/**

## Results

![Concentration in the domain](diffusion_1D_transient/c.png)(width="800" height="600")

~~~gnuplot Concentration profile at t=60s
set xlabel 'x [m]'
set ylabel 'c [mol]'
plot 'out' i 5:5 u 1:4 w l t 'exact lvl10',\
'out' i 1:1 u 1:3 w p t 'Basilisk lvl 8',\
'out' i 3:3 u 1:3 w p t 'Basilisk lvl9',\
'out' i 5:5 u 1:3 w p t 'Basilisk lvl10'
~~~

~~~gnuplot Concentration profile at t=1 and 60s
set xlabel 'x [m]'
set ylabel 'c [mol]'
plot 'out' i 4:4 u 1:4 w l t 'exact lvl10 1s',\
'out' i 5:5 u 1:4 w l t 'exact lvl10 60s',\
'out' i 4:4 u 1:3 w p t 'Basilisk lvl10 1s',\
'out' i 5:5 u 1:3 w p t 'Basilisk lvl10 60s'
~~~
*/
