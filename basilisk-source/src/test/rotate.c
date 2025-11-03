/**
# Rotation of a circular interface

~~~gnuplot VOF facets for N = 64
set size ratio -1
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
plot 'interface' w l t ''
~~~

~~~gnuplot Error field for N = 64
reset
set term @PNG
set output "error.png"
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,\
     0.875 0.9333 0 0, 1 0.498 0 0 )
set size ratio -1
splot [0.1:0.4][-0.15:0.15]'out' t ""
~~~

~~~gnuplot Error convergence
reset

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

fc(x)=ac+bc*x
fit fc(x) 'clog' u (log($1)):(log($4)) via ac,bc
fc2(x)=ac2+bc2*x
fit fc2(x) 'clog' u (log($1)):(log($2)) via ac2,bc2

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]
plot 'log' u 1:4 t 'max (adaptive)', exp(f(log(x))) t ftitle(a,b), \
     'clog' u 1:4 t 'max (constant)', exp(fc(log(x))) t ftitle(ac,bc), \
     'log' u 1:2 t 'norm1 (adaptive)', exp(f2(log(x))) t ftitle(a2,b2), \
     'clog' u 1:2 t 'norm1 (constant)', exp(fc2(log(x))) t ftitle(ac2,bc2)
~~~
*/

#include "advection.h"
#include "vof.h"

scalar c[];
scalar * interfaces = {c}, * tracers = NULL;
int MAXLEVEL;

int main()
{
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  for (MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run ();
  }
}

#define circle(x,y) (sq(0.1) - (sq(x-0.25) + sq(y)))

event init (i = 0)
{
  fraction (c, circle(x,y));
}

#define end 0.785398

event velocity (i++) {
#if TREE
  double cmax = 1e-2;
  adapt_wavelet ({c}, &cmax, MAXLEVEL, list = {c});
#endif

  double a = -8.;
  trash ({u});
  foreach_face(x) u.x[] = - a*y;
  foreach_face(y) u.y[] =   a*x;
}

event logfile (t = {0,end}) {
  stats s = statsf (c);
  fprintf (stderr, "# %f %.12f %f %g\n", t, s.sum, s.min, s.max);
}

event interface (t += end/10.) {
  if (N == 64) {
    static FILE * fp = fopen ("interface", "w");
    output_facets (c, fp);
  }
}

event field (t = end) {
  scalar e[];
  fraction (e, circle(x,y));
  foreach()
    e[] -= c[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  if (N == 64)
    output_field ({e}, stdout, N, linear = false);
}
