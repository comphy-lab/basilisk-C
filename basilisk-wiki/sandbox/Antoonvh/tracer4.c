/**
# Tracer advection in a accelerating frame of reference

![The concergence test](tracer4/s.mp4)

![Variance test (advection of noise)](tracer4/n.mp4)

~~~gnuplot Convergence
set xr [8:512]
set yr [1e-5:10]
set grid
set logscale x 2
set logscale y
set size square
set xlabel 'N'
set ylabel 'L_1'
plot 'out' t 'Data', 1e6*x**(-4) t '4th order' 
~~~

~~~gnuplot Advection of noise
reset
set xr [8:512]
set yr [0:40]
set grid
set logscale x 2
set size square
set xlabel 'N'
set ylabel 'Sum of squares'
plot 'out' u  1:3 t 'Data', 33.333 t 'Perfect noise' 
~~~
 */
#include "nsf4t.h"
#include "view.h"

scalar s[], n[], * tracers = {s, n};
double ax = 0.1, ay, uo = 0.5, tend;

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 10;
  X0 = Y0 = -L0/2;
  tend = (-uo + sqrt(sq(uo) - 2*ax*-L0))/(ax); // ABC formula
  ay = 2*L0/sq(tend);
  const vector av[] = {ax, ay};
  a = av;
  for (N = 16; N <= 256; N *= 2)
    run();
}

event init (t = 0) {
  foreach_face(x)
    u.x[] = uo;
  foreach_vert() {
    s[] = exp(-sq(x) - sq(y));
    n[] = noise();
  }
}

event mover (t += 0.5) {
  squares ("s", min = -0.1, max = 1.1);
  save ("s.mp4");
  squares ("n", min = -1, max = 1, map = cool_warm);
  save ("n.mp4");
  // There is an issue with these:
  // output_ppm (s, file = "s.png", n = 300);
  // output_ppm (n, file = "n.mp4", n = 300, min = -1, max = 1);
}

event stop (t = tend) {
  event ("mover");
  double e = 0, E = 0;
  foreach_vert() {
    e += dv()*fabs(s[] - exp(-sq(x) - sq(y)));
    E += dv()*sq(n[]);
  }
  printf ("%d %g %g\n", N, e, E);
}


