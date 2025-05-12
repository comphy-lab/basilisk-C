/**
# 2D Heat equation
we solve T(x, y, t)
$$ \frac{\partial T}{\partial t} = c^{2} (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2})$$

The discretization form of the equation:
$$\frac{\partial^2 T}{\partial x^2} |_{(i, j)}  = \frac{1}{\Delta_x^{2}} [T_{i+1, j} - 2 T_{i, j} + T_{i-1, j}] $$
$$\frac{\partial^2 T}{\partial y^2} |_{(i, j)}  = \frac{1}{\Delta_y^{2}} [T_{i, j+1} - 2 T_{i, j} + T_{i, j-1}] $$
*/

#include "run.h"
#include "diffusion.h"
#include "view.h"

scalar T1[], dTx[], dTy[];
scalar T2[];
scalar T3[];
double dt;
double alpha;
#define EPS 1.1

int main()
{
  L0 = 20;
  X0 = -L0 / 2.;
  Y0 = -L0 / 2.;
  N = 1 << 6;
  DT = sq(L0 / N) / 4;
  run();
}

T3[left] = dirichlet(3.);
T3[top] = dirichlet(0.);

event init(t = 0)
{
  foreach ()
  {
    T1[] = 10. * (sqrt(sq(x) + sq(y)) < EPS) / 2;
    T2[] = y >= 0 ? 3 : 0;
    T3[] = y >= 0 ? 3 : 0;
  }
  boundary({T1, T2, T3});
}

event integration(i++)
{
  double alpha = 1;                         // uniform material
  DT = sq(L0 / (1 << grid->maxdepth)) / 4.; // smaller time step when the grid is refined
  double dt = DT;
  dt = dtnext(dt);
  foreach ()
  {
    dTx[] = (T2[1, 0] - 2 * T2[0, 0] + T2[-1, 0]) / sq(Delta);
    dTy[] = (T2[0, 1] - 2 * T2[0, 0] + T2[0, -1]) / sq(Delta);
  }
  foreach ()
    T2[] += dt * alpha * (dTx[] + dTy[]);
  boundary({T2});
}

mgstats mgb1, mgb3;
event Diffusion(i++)
{
  DT = sq(L0 / (1 << grid->maxdepth)) / 3.;
  double dt = DT;
  mgb1 = diffusion(T1, dt);
  mgb3 = diffusion(T3, dt);
  // boundary({T1, T3});
}

/**
Let try different output options:

* directly print
* output_field
* output_ppm
* save movies

*/

event printdata(t = 0; t <= 1; t += 0.05)
{
  static FILE *fp = fopen("output_H2DT1", "w");
  for (double y = -L0 / 2; y < L0 / 2; y += L0 / N)
  {
    fprintf(fp, "%g %g\n", y, interpolate(T1, 0, y));
  }
  fprintf(fp, "\n \n");
  fflush(fp);

  static FILE *fp2 = fopen("output_H2DT2", "w");
  for (double y = -L0 / 2; y < L0 / 2; y += L0 / N)
  {
    fprintf(fp2, "%g %g\n", y, interpolate(T2, 0, y));
  }
  fprintf(fp2, "\n \n");
  fflush(fp2);

  static FILE *fp3 = fopen("output_H2DT3", "w");
  for (double y = -L0 / 2; y < L0 / 2; y += L0 / N)
  {
    fprintf(fp3, "%g %g\n", y, interpolate(T3, 0, y));
  }
  fprintf(fp3, "\n \n");
  fflush(fp3);
}

event outmatrix(t = 0.2)
{
  FILE *fp1 = fopen("outputfield_H2DT1", "w");
  output_field({T1}, fp1);

  FILE *fp2 = fopen("outputfield_H2DT2", "w");
  output_field({T2}, fp2);

  FILE *fp3 = fopen("outputfield_H2DT3", "w");
  output_field({T3}, fp3);
}

event movie(t = 0; t <= 0.4; t += 0.01)
{
  output_ppm(T1, file = "T1_ap.mp4", min = 0, max = 5, linear = true, map = cool_warm);
  output_ppm(T2, file = "T2_ap.mp4", min = 0, max = 5, linear = true, map = cool_warm);
  output_ppm(T3, file = "T3_ap.mp4", min = 0, max = 5, linear = true, map = cool_warm);
}

event mov(i++)
{
  squares("T1", linear = true, min = -0.1, max = 5, map = cool_warm);  // not sure why my cells overlay on Tsquare
  cells();
  save("grid.mp4");
}

event adapt(i++)
{
  adapt_wavelet({T1, T2, T3}, (double[]){3e-1, 3e-1, 3e-1}, minlevel = 4, maxlevel = 8);
}


/**
let's check some profiles
~~~gnuplot
file="H2DT"

set terminal png size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set multiplot layout 1,3
set xlabel "x"
set ylabel "T"
p[][] 'output_'.file.'1' u ($1):($2) t 'IC1 BC1' w l
p[][] 'output_'.file.'2' u ($1):($2) t 'IC2 BC2' w l
p[][] 'output_'.file.'3' u ($1):($2) t 'IC2 BC2' w l
unset multiplot
~~~
*/



/**
this is funny to see
<table style="margin-left:auto;margin-right:auto;width:100%;">
<tr><td>
![T1](heat2D/T1_ap.mp4)(width="100%")
</td><td>
![T2](heat2D/T2_ap.mp4)(width="100%")
</td><td>
![T3](heat2D/T3_ap.mp4)(width="100%")
</td>
</tr>
</table>
*/

/**
let's check the grid
![grid](heat2D/grid.mp4)(width="60%")
So the grid changes for all variables not just one of them!
*/
