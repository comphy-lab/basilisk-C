/**
# Ekman flow in 2D

We test the `force_geo.h` header file in 2D
 */

#include "navier-stokes/centered.h"
#define U_GEO 1.
#include "force_geo.h"
scalar * tracers;      //We need to declare the tracers

uz[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);

FILE * fp;

int main() {
  fp = fopen ("vy", "w");
  periodic (left);
  f_cor = 1.;
  mu = unityf;
  L0 = 25.;
  for (TOLERANCE = 1e-2; TOLERANCE >= 1e-6; TOLERANCE /= 10.)
    run();
}

event init (t = 0) {
  foreach()
    u.x[] = U_GEO;
}

event adapt (i++)
  adapt_wavelet ({u.x, uz}, (double[]){0.005, 0.005}, 8);

#include "profile5c.h"
scalar u_test[];
event end (t < HUGE; i += 10) {
  if (change (u.x, u_test) < 2e-3) {
    fprintf (fp, "%g %g\n", TOLERANCE, normf(u.y).avg);
    if (fabs (TOLERANCE - 1e-3) < 1e-4)
      profile ({u.x, uz}, fname = "prof");
    return 1;
  }
}

/**
## Results

The approximate solution (slowly) converges towards the analytical
Ekman spiral over time.

~~~gnuplot profiles
set xr [0 : 25]
set xlabel 'y'
set ylabel 'velocities'
  plot 'prof' u 1:2 t 'u.x', 'prof' u 1:3 t 'uz',\
          1 - exp(-0.707*x)*cos(0.707*x), exp(-0.707*x)*sin(0.707*x)
~~~

The effect of the `TOLERANCE` on the average absolute vertical velocity

~~~gnuplot Vertical velocity
reset
set xr [0.0000001:0.1]
set xlabel 'TOLERANCE'
set ylabel 'average  abs(u.y)'
set key off
set logscale xy
plot 'vy' u 1:2
~~~
 */
