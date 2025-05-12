/**
# Taylor-Green vortices

The [`nsrk.h`](nsrk.h) solver is more accurate than this [Centered
solver](/src/navier-stokes/centered.h) and the [all-Mach
solver](/src/all-mach.h) for [this test](/src/test/taylor-green.c).

~~~gnuplot Error Norms
set xlabel 'Spatial resolution'
set ylabel 'Error norms'
set cbrange [1:1]
set logscale
set grid
set xtics 16,2,256
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit [4:] f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit [4:] fm(x) 'log' u (log($1)):(log($4)) via am,bm
set xrange [16:512]
set pointsize 1
plot exp (f2(log(x))) t ftitle(a2,b2), \
     exp (fm(log(x))) t ftitle(am,bm),  \
     'log' u 1:3 t '|e|_2', \
     'log' u 1:4 t '|e|_{max}' lc 0
~~~

~~~gnuplot The Divergence is excellent because `it_project=1`
reset
set xlabel 'Time'
set ylabel 'Maximum divergence'
set cbrange [1:1]
set grid
set logscale y
set xrange [0:2]
set yrange [1e-8:]
plot '< grep "^32 " out' u 3:4 w l t '32^2', \
     '< grep "^64 " out' u 3:4 w l t '64^2', \
     '< grep "^128 " out' u 3:4 w l t '128^2', \
     '< grep "^256 " out' u 3:4 w l t '256^2'
~~~

~~~gnuplot Equivalent Reynolds number
reset
set xlabel 'Resolution'
set ylabel 'Equivalent Reynolds number'
set grid
set cbrange [1:1]
set logscale
set xtics 16,2,256
set xrange [16:512]
f(x)=a*exp(-b*x)
Re(b)=1./(b/(4.*(2.*pi)**2))

set print "Re"
fit [1:] f(x) '< grep "^32 " out' u 3:5 via a,b
print 32,Re(b)
fit [1:] f(x) '< grep "^64 " out' u 3:5 via a,b
print 64,Re(b)
fit [1:] f(x) '< grep "^128 " out' u 3:5 via a,b
print 128,Re(b)
fit [1:] f(x) '< grep "^256 " out' u 3:5 via a,b
print 256,Re(b)

plot 'Re' w lp t 'nsrk'
~~~

 */
#include "nsrk.h"

int main() {
  origin (-0.5,-0.5);
  foreach_dimension()
    periodic (right);
  it_project = 1;
  for (N = 32; N <= 256; N *= 2)
    run();
}

event init (t = 0) {
  CFL = 0.9;
  foreach() {
    u.x[] = - cos(2.*pi*x)*sin(2.*pi*y);
    u.y[] =   sin(2.*pi*x)*cos(2.*pi*y);
    p[]   = - (cos(4.*pi*x) + cos(4.*pi*y))/4.;
  }
}

event logfile (i++) {
   scalar div[], ke[];
  foreach() {
    div[] = (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
    ke[] = sq(u.x[]) + sq(u.y[]);
  }
  printf ("%d %d %g %g %g\n", N, i, t, normf(div).max, statsf(ke).sum);
}

event error (t = 2) {
  scalar e[];
  foreach() {
    double u0 = - cos(2.*pi*x)*sin(2.*pi*y);
    double v0 =   sin(2.*pi*x)*cos(2.*pi*y);
    e[] = norm(u) - sqrt(sq(u0) + sq(v0));
  }
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
