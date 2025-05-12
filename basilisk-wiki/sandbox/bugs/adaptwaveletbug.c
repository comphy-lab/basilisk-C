/**
This case includes the axisymmetric solver combining with Navier-Stokes, two-phase flow and surface tension as well
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

/**
Introduce some constants of the case,
*/

int LEVELmax = 9, LEVELmin = 5;

/**
The program starts here,
*/

int main()
{
  size(10.0);
  init_grid(16);
  rho1 = 1000.0;
  rho2 = 1.0;
  mu1 = 5.0e-4;
  mu2 = 2.0e-5;
  f.sigma = 1.6e-2;
  run();
  return 1;
}

/**
As it can be seen, the initial condition is showing a drop with 10% gap for the refinement. The drop velocity is also set here.
*/

event init(i = 0)
{
  double x0 = 7.50;
  refine(sq(x - x0) + sq(y) + sq(z) < sq(1.1) && sq(x - x0) + sq(y) + sq(z) > sq(0.9) && level < LEVELmax);
  foreach()
  {
    if (sq(x - x0) + sq(y) + sq(z) < sq(1.0))
    {
      f[] = 1.0;
      u.x[] = -10.0;
    }
    else
      f[] = 0.0;
  }
}

/**
The adapt_wavelet is used to construct a refined mesh close to the interface by choosing the "f" as the variable and "0.001" as the tolerance.
*/

event adapt(i++)
{
  adapt_wavelet({f}, (double[]) {0.001}, maxlevel = LEVELmax, minlevel = LEVELmin);
}

/**
And, finally, for the iteration showing and GFSView output we will use the following "events".
*/

event showiteration(i += 50)
{
  printf("i[%06d], dt[%e], t[%.2f]\r\n", i, dt, t);
}

event gfsview(t += 0.1)
{
  static FILE * fp;
  char name[500];
  sprintf(name, "GFS[%05d]-t[%.2f].gfs", i, t);
  fp = fopen(name, "w");
  output_gfs(fp, translate = true);
  printf("write data file!\r\n");
  fclose(fp);
}

/**
The final time can be set here as well.
*/

event end(t = 0.50)
{
  printf("i[%06d], dt[%e], t[%.2f]\r\nEND OF RUN!\r\n", i, dt, t);
}

/**
The results has to show a drop moving with constant velocity; therefore, the adapt_wavelet has to move with the interface accordingly. Nevertheless, the refinement is not following the interface in some cells.
*/