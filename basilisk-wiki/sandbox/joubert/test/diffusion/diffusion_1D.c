/**
# 1D Diffusion test steady state with imposed concentration boundary condition. */

//Geometric parameters of the drop
#define H (0.27)

// Physical parameters
#define Sctw (1)
#define NUw (1e-6)

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

#define MAXTIME (36450)
#define cexact(x) (1+((x)*(-1))/0.27)

scalar c[], cold[];

c[left] = dirichlet(1);
c[right]  = dirichlet(0);
  
int maxlevel;
int minlevel = 6;
double dt,errmax;
mgstats mgd1;

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

event time_integration (i++) {
 dt = dtnext (DT);
  const face vector D[] = {NUw/Sctw};
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


event end (t = MAXTIME) {
  scalar dce[];
  foreach()
    dce[] = c[]-cold[];
  boundary({dce});
  output_ppm (c, file = "c.png",spread=2, linear = true, n=512);
  output_ppm (dce, file = "difference_c.png",spread=2, linear = true, n=512);
  foreach()
    fprintf (fout, "%g %d %g %g\n", x, N, c[], cexact(x));
    fflush (fout);
}

/**

## Results

![Concentration on the domain](diffusion_1D/c.png)(width="800" height="600")

~~~gnuplot Evolution of concentration along the domain
set xlabel 'x [m]'
set ylabel 'c [mol]'
plot 'out' i 2:2 u 1:4 w l t 'exact lvl10',\
'out' i 0:0 u 1:3 w p t 'Basilisk lvl 8',\
'out' i 1:1 u 1:3 w p t 'Basilisk lvl9',\
'out' i 2:2 u 1:3 w p t 'Basilisk lvl10'
~~~

*/
