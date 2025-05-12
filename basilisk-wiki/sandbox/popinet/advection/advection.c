/**
# Advection of a scalar field

Using either the second-order Bell--Colella--Glaz scheme or the
5th-order WENO scheme. */

#include "grid/multigrid.h"
#include "advection.h"
#include "view.h"

static double func (double x, double y, double offsetx, double offsety) {
  if (sq(x - offsetx) <= 0.0225 && sq(y - offsety) <= 0.0225)
    return (pow(1. - sq((x - offsetx)/0.15), 7)*
	    pow(1. - sq((y - offsety)/0.15), 7));
  else
    return 0;
}

static double integral (double x, double y, double Delta,
			double offsetx, double offsety)
{
  double sum = 0.;
  double q = Delta/2.*sqrt(3./5.);
  sum += 5.*(5.*func (x - q, y - q, offsetx, offsety) +
	     8.*func (x, y - q, offsetx, offsety) +
	     5.*func (x + q, y - q, offsetx, offsety));
  sum += 8.*(5.*func (x - q, y, offsetx, offsety) +
	     8.*func (x, y, offsetx, offsety) +
	     5.*func (x + q, y, offsetx, offsety));
  sum += 5.*(5.*func (x - q, y + q, offsetx, offsety) +
	     8.*func (x, y + q, offsetx, offsety) +
	     5.*func (x + q, y + q, offsetx, offsety));
  return sum/sq(18.);
}

void face_velocity (face vector u, double t)
{
  foreach_face(x)
    u.x[] =  1.5/(pi*Delta)*sin(2.*pi*t/5.)*cos(pi*x)*
    (cos(pi*(y + Delta/2.)) - cos(pi*(y - Delta/2.)));
  foreach_face(y)
    u.y[] = - 1.5/(pi*Delta)*sin(2.*pi*t/5.)*cos(pi*y)*
    (cos(pi*(x + Delta/2.)) - cos(pi*(x - Delta/2.))); 
  boundary((scalar *){u});
}

#if !WENO
event velocity (i++) {
  face_velocity (u, t);
}
#endif

scalar f[], * tracers = {f};

int main()
{
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
  for (N = 64; N <= 512; N *= 2)
    run();
}

event init (i = 0)
{ 
  foreach()
    f[] = integral (x, y, Delta, - 0.2, - 0.236338);
  boundary({f});
}

event logfile (t = {0,5}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %ld %.12f %g %g\n", t, grid->tn, s.sum, s.min, s.max);
}

event deformed (t = 2.5) {
  view (fov = 19.35, width = 400, height = 400);
  clear();
  squares ("f", linear = true);
  save ("f.png");
}

event field (t = 5) {
  scalar e[];
  foreach()
    e[] = f[] - integral (x, y, Delta, -0.2, -0.236338);
  boundary ({e});
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);

  view (fov = 19.35, width = 400, height = 400);
  clear();
  squares ("e", linear = true);
  save ("e.png");
}

/**
## Results

![Tracer field at $t=2.5$](advection/f.png)

![Error field at $t=5$ (BCG)](advection/e.png)

![Error field at $t=5$ (WENO)](advection-weno/e.png)

~~~gnuplot
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit [4.5:] f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit [4.5:] f2(x) '' u (log($1)):(log($2)) via a2,b2
f3(x)=a3+b3*x
fit [4.5:]f3(x) '../advection-weno/log' u (log($1)):(log($4)) via a3,b3
f4(x)=a4+b4*x
fit [4.5:]f4(x) '' u (log($1)):(log($2)) via a4,b4
set xlabel 'Resolution'
set ylabel 'Error'
set logscale
set xrange [32:1024]
set xtics 32,2,1024
set grid ytics
set cbrange [1:2]
set ytics format "%.0e"
set key bottom left
plot 'log' u 1:4 t 'max (BCG)', exp(f(log(x))) t ftitle(a,b),      \
     '' u 1:2 t 'norm1 (BCG)', exp(f2(log(x))) t ftitle(a2,b2),    \
     '../advection-weno/log' u 1:4 t 'max (WENO)',                 \
     exp(f3(log(x))) t ftitle(a3,b3),			           \
     '' u 1:2 t 'norm1 (WENO)', exp(f4(log(x))) t ftitle(a4,b4)
~~~

~~~gnuplot Error
reset
set xlabel 'Computing time (sec)'
set ylabel 'Error'
set logscale
set grid ytics
set ytics format "%.0e"
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $2,$4}' < log > error
! awk '($1 == "#") {print $5}' < ../advection-weno/out > time-weno
! awk '($1 != "#") {print $2,$4}' < ../advection-weno/log > error-weno
set key bottom left
plot '< paste time error' u 1:3 w lp t 'max (BCG)',            \
     '< paste time-weno error-weno' u 1:3 w lp t 'max (WENO)', \
     '< paste time error' u 1:2 w lp t 'norm1 (BCG)',	       \
     '< paste time-weno error-weno' u 1:2 w lp t 'norm1 (WENO)'
~~~
*/
