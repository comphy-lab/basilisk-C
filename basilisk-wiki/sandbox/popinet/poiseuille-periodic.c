/**
# Poiseuille flow of a yield-stress fluid
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "yield.h"

double tau_y = 0.2;

int main() {
  periodic (top);
  size (0.5);
  DT = HUGE [0];
  
  stokes = true;
  TOLERANCE = 1e-5;
  
  u.t[left] = dirichlet(0);

  const face vector g[] = {0.,1.};
  a = g;
  const face vector muc[] = {1.,1.};
  mu = muc;
  const face vector tauc[] = {tau_y,tau_y};
  tau_0 = tauc;

  for (N = 8; N <= 64; N *= 2)
    run();
}

double analytical (double x)
{
  double h = 0.5, dpdy = 1., mu = 1.;
  double U = sq(tau_y - h*dpdy)/(2.*mu*dpdy);
  double X = h - tau_y/dpdy;
  return (x < X ? dpdy*(h*x - sq(x)/2.) - tau_y*x : U);
}

scalar un[];

event logfile (t += 0.1; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-7)
    return 1; /* stop */
  //  fprintf (stderr, "%g %d %g %d %d\n", t, i, du, mgp.i, mgu.i);
}

event profile (t = end) {
  if (N == 16) {
    foreach()
      fprintf (stdout, "%g %g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[],
	       x - Delta/2., (u.y[] - u.y[-1])/Delta);
    printf ("\n");
  }
  scalar e[];
  foreach()
    e[] = u.y[] - analytical(x);
  norm n = normf (e);
  fprintf (ferr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

/**
~~~gnuplot Velocity profile: comparison between numerical and analytical solutions
set xlabel 'x'
set ylabel 'u.y'
h = 0.5
dpdy = 1.
tau_y = 0.2
mu = 1.
U = (tau_y - h*dpdy)**2/(2.*mu*dpdy)
X = h - tau_y/dpdy
u(x) = x < X ? dpdy*(h*x - x**2/2.) - tau_y*x : U
set xrange [0:0.5]
plot 'out' u 1:4 w p t 'basilisk', u(x) w l t 'analytical'
~~~

~~~gnuplot Shear stress: comparison between numerical and analytical solutions
set xlabel 'x'
set ylabel 'du.y/dx'
plot 'out' u 6:7 w p t 'basilisk', (x<X?X-x:0) w l t 'analytical'
~~~

~~~gnuplot Convergence of the error on u.y
reset
set xlabel 'Resolution'
set ylabel 'Error norms'
set key bottom left
set logscale
set cbrange [1:2]
set xtics 16,2,128
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' u (log($1)):(log($4)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), \
     exp (f2(log(x))) t ftitle(a2,b2), \
     'log' u 1:4 t '|h|_{max}' ps 1.5, \
     'log' u 1:3 t '|h|_2' ps 1.5
~~~

## See also

* [Free surface flow of a Bingham
   fluid](/sandbox/M1EMN/Exemples/bingham_simple.c)
*/
