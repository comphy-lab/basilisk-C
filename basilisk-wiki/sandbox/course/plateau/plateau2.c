/**
# Rayleigh-Plateau instability

* vary epsilon from 0.1 -> 0.001 and see if the growth rate changes
* vary R and see how the growth rate changes (0.4 stable, 0.3 unstable)
* does the spatial resolution change things?
* influence of density ratio?
* points on the growth rate curve

![](plateau2/f.mp4)
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double R = 0.2, epsilon = 0.1;

int main() {
  //  rho2 = 0.01;
  f.sigma = 1.;
  run();
}

event init (t = 0) {
  double k = pi;
  fraction (f, R*(1. + epsilon*cos(k*x)) - y);
}
#if 0
event gfsview (i += 10; t <= 1) {
  static FILE * fp = popen ("gfsview2D -s plateau2.gfv", "w");
  output_gfs (fp);
}
#endif

event mov (t += 0.1; t<= 1) {
  output_ppm (f, file = "f.mp4", opt = "-r 5", n = 300);
}

event logfile (i++) {
  scalar hy[];
  foreach() {
    if (f[] > 1e-3 && f[] < 1 - 1e-3) {
      coord m = mycs (point, f), fc;
      double alpha = plane_alpha (f[], m);
      plane_area_center (m, alpha, &fc);
      hy[] = y + Delta*fc.y;
    }
    else
      hy[] = nodata;
  }
  stats s = statsf (hy);
  fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

/**
~~~gnuplot Evolution of maximum radius with time
set xlabel 'Time'
set ylabel 'Rmax'
set logscale y
set grid
R=0.2
epsilon=0.1

lambda=2.
k=2.*pi/lambda
x=k*R
sigma=1
rho=1
omega0=sqrt(sigma/(2.*rho*R**3))

# x = 0.63, 3.4/omega0 = 0.43 instead of 0.544

omega(x)=sqrt(x**2*(1-x**3))

plot [0:0.8]'log' u 1:($3 - R)/(R*epsilon) w l, exp(3.4*x)
~~~
*/
