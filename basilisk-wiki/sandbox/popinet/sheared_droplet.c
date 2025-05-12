/**
# Surface transport as "shell transport"

![A droplet with a "shell" in a Couette flow.](sheared_droplet/movie.mp4)

~~~gnuplot "Surface concentration" at $t = 0$ and $t = 1.5$.
unset key
set xlabel 'Angle'
set xtics ("-π" -pi, "-π/2" -pi/2, "0" 0, "π/2" pi/2, "π" pi)
set ylabel 'Concentration'
plot [-pi:pi] 'log' u 1:2 w p
~~~
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#define Ca 0.3     // Capillary number
#define Re 0.3     // Reynold number
#define We (Ca*Re) // Weber number
#define MUr 0.5    // ratio of outer(matrix) to inner(drop) viscosity
#define M 1.       // ratio of outer to inner density

int MAXLEVEL = 8;
int MINLEVEL = 6;

scalar f1[], f2[]; // these are just for visualisation
scalar f3[]; // this is used for "surface transport"

u.t[top] = dirichlet (y);
u.t[bottom] = dirichlet (y);

int main() {
  L0 = 16.;
  origin (-L0/2, -L0/4.);
  periodic (right);
  DT = 1e-1;

  mu1 = MUr/Re;
  mu2 = 1./Re;
  rho1 = M;
  rho2 = 1.;
  f.sigma = 1./We;

  init_grid (1 << MAXLEVEL);
  interfaces = list_concat (NULL, {f, f1, f2, f3});
  run();
}

/**
When epsilon is small enough, the results are independent from its value. */

const double epsilon = 1e-6;

event init (i = 0) 
{
  mask (y > 4 ? top : none);

  fraction (f, 1. - (sq(x) + sq(y)));
  fraction (f1, sq(1. + 0.02) - (sq(x) + sq(y)));
  fraction (f2, sq(1. - 0.02) - (sq(x) + sq(y)));
  fraction (f3, sq(1. + epsilon) - (sq(x) + sq(y)));
  foreach()
    u.x[] = y;
}

event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);
}

event movie (i += 2)
{
  view (quat = {0.000, 0.000, 0.000, 1.000},
	fov = 30, near = 0.01, far = 1000,
	tx = -0.003, ty = -0.005, tz = -0.340,
	width = 800, height = 600);
  vectors (u = "u", scale = 0.01);
  draw_vof (c = "f", lw = 2);
  draw_vof (c = "f1", lc = {1,0,0}, lw = 2);
  draw_vof (c = "f2", lc = {0,0,1}, lw = 2);
  save ("movie.mp4");
}

event profiles (t = {0, 1.5})
{
  foreach()
    if (f[] > 1e-6 && f[] < 1. - 1e-6 &&
	f3[] > 1e-6 && f3[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      double alpha3 = plane_alpha (f3[], n);
      double d = Delta*(alpha3 - alpha)/sqrt(n.x*n.x + n.y*n.y);
      
      plane_area_center (n, alpha, &p);
      double theta = atan2 (y + Delta*p.y, x + Delta*p.x);
      fprintf (stderr, "%g %g\n", theta, d/epsilon);
    }
}

event cleanup (t = end)
{
  free (interfaces);
}
