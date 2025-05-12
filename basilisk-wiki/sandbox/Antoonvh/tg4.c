/**
# Taylor-Green vortices

Convergence test for an inviscid and viscous/diffusive medium.

![Advecting vortices (vorticity)](tg4/o.mp4)

![Advecting scalar](tg4/s.mp4)


~~~gnuplot Error Velocity
set xr [5: 100]
set xlabel 'Spatial resolution (N)'
set ylabel 'Error'
set size square
set grid
set logscale x 2
set logscale y
plot 'out' u 1:2 t 'L_1 inviscid' ps 2, '' u 1:3 t 'L_{inf} inviscid' ps 2, 1e3*x**(-4) t '4th order', \
'log' u 1:2 t 'L_1 viscous' ps 2, '' u 1:3 t 'L_{inf} viscous'
~~~

~~~gnuplot Error tracer
reset
set xr [5: 100]
set xlabel 'Spatial resolution (N)'
set ylabel 'Error'
set size square
set grid
set logscale x 2
set logscale y
plot 'out' u 1:5 t 'L_1 inviscid' ps 2, '' u 1:6 t 'L_{inf} inviscid' ps 2, 1e3*x**(-4) t '4th order', \
'log' u 1:5 t 'L_1 diffusive' ps 2, '' u 1:6 t 'L_{inf} diffusive'
~~~
 */
#include "nsf4t.h"
scalar s[], * tracers = {s};

double muv = 0;

int main() {
  X0 = Y0 = -L0/2;
  foreach_dimension()
    periodic (right);
  for (N = 8; N <= 64; N *= 2)
    run();
  muv = 0.005;
  const scalar muc[] = muv;
  nu = kappa = muc;
  for (N = 8; N <= 64; N *= 2)
    run();
  
}

double u_x (double x, double y) {
  return -cos(2.*pi*x)*sin(2.*pi*y)*exp(-2.*muv*sq(2.*pi)*t) + 1.;
}

double u_y (double x, double y) {
  return  sin(2.*pi*x)*cos(2.*pi*y)*exp(-2.*muv*sq(2.*pi)*t) + 0.5;
}

double s_a (double x, double y) {
  return  cos(2.*pi*x)*cos(2.*pi*y)*exp(-2*muv*sq(2.*pi)*t);
}

event init (t = 0) {
  TOLERANCE = 1e-5;
  foreach_face() 
    u.x[] = Gauss6_x (x, y, Delta, u_x);
  foreach_vert()
    s[] = s_a(x, y);
}

event mov (i += 5) {
  scalar omg[];
  vorticityf (u, omg);
  output_ppm (omg, file = "o.mp4", n = 300, min = -10, max = 10);
  output_ppm (s, file = "s.mp4", n = 300,
	      map = cool_warm, min = -1, max = 1);
}

event errors (t = 2) {
  if (i % 5 != 0)
    event ("mov");
  double e = 0, es = 0, E = 0;
  double em = -1, esm = -1;
  foreach_face() {
    E  += sq(Delta)*sq(u.x[]);
    double el = fabs(Gauss6_x (x, y, Delta, u_x) - u.x[]);
    e  += sq(Delta)*el;
    if (el > em)
      em = el;
  }
  foreach_vert() {
    double el = fabs(s_a(x, y) - s[]);
    es += sq(Delta)*el;
    if (el > esm)
      esm = el;
  }
  fprintf (muv ? stderr : stdout, "%g %g %g %g %g %g\n", sqrt(grid->tn), e, em, E, es, esm);
  return 1;
}
