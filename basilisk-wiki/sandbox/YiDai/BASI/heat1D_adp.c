/**
# solver 1D heat equation. 
$$\partial_{t}T = \partial_{xx}T$$

Here we use adaptive grid functionality 

*/

#include "grid/bitree.h"
#include "diffusion.h"
#include "run.h"


// two initial condition
scalar T1[], dT1[];
scalar T2[];
double dt;
#define EPS 1

int main()
{
  L0 = 20;
  X0 = -L0 / 2.;
  N = 1 << 4;
  DT = sq(L0 / N) / 3; // smaller than von neumann stability
  run();
}

event init(t = 0)
{
  foreach ()
  {
    T1[] = 1. / EPS * (fabs(x) < EPS) / 2;
    T2[] = x >= 0 ? 3 : 0;
  }
  boundary({T1, T2});
}

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata(t = 0; t <= 1; t += 0.2)
{
  static FILE *fp = fopen("output_H1DA", "w");
  foreach ()
    fprintf(fp, "%g %g %g %g %g\n", x, T1[], T2[], dT1[], t);
  fprintf(fp, "\n\n");
  fflush(fp);
}

event integration(i++)
{
  // DT = sq(L0 / (1 << grid->maxdepth)) / 4.; // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext(dt);
  foreach ()
  {
    dT1[] = (T1[1, 0] - 2 * T1[0, 0] + T1[-1, 0]) / sq(Delta);
  }
  foreach ()
  {
    T1[] += dt * dT1[];
  }
  boundary({T1});
}

mgstats mgb;
event Diffusion(i++)
{
  DT = sq(L0 / (1 << grid->maxdepth)) / 3.;
  double dt = DT;
  mgb = diffusion(T2, dt);
  boundary({T2});
}

event adapt(i++)
{
  adapt_wavelet({T1, T2}, (double[]){0.5}, minlevel = 2, maxlevel = 6);
}

/**
~~~gnuplot
set terminal png size 800,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

file="H1DA"
set multiplot layout 1,2
set xlabel "x"
set ylabel "T"
p[][] 'output_'.file u ($1):($2) t 'IC1' w lp
p[][-0.5:3.5] 'output_'.file u ($1):($3) t 'IC2' w lp
unset multiplot
~~~
*/

/**
let check out the adaptive grid at t = 0.6. ha still hard to see
~~~gnuplot
set terminal png size 800,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

file="H1DA"
set multiplot layout 1,2
set xlabel "x"
set ylabel "T"
p[][] 'output_'.file u ($1):($5==0.6?$2:1/0) t 'IC1' w lp
p[][-0.5:3.5] 'output_'.file u ($1):($5==0.6?$3:1/0) t 'IC2' w lp pointtype 21
unset multiplot
~~~
*/