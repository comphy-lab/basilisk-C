/**
# Taylor--Green vortices

[Taylor--Green
vortices](http://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex)
are one of the few exact non-trivial solutions of the incompressible
Euler equations. In this test case, we use this solution as initial
condition and check whether the numerical scheme can respect the
balance between non-linear advection terms and pressure
gradients. Numerical diffusion will in particular introduce
dissipation. This dissipation can be quantified and is a useful
measure of the accuracy of the numerical scheme.

We solve the incompressible Euler equations on a Cartesian (multi)grid
either with the centered Navier-Stokes solver or with the "all Mach"
solver (in incompressible mode). */

#include "grid/multigrid.h"
#if ALL_MACH
# include "all-mach.h"
# include "bcg.h"
# define u q

event tracer_advection (i++)
  advection ((scalar *){q}, uf, dt, (scalar *){g});
#else
# include "navier-stokes/centered.h"
#endif

int main() {

  /**
  Space and time are dimensionless. The domain is unity, centered on
  the origin and periodic in all directions. */

  size (1. [0]);
  DT = HUGE [0];
  origin (-0.5,-0.5);
  foreach_dimension()
    periodic (right);

  /**
  We check convergence with spatial resolution from 32^2^ to
  256^2^. */
  
  for (N = 32; N <= 256; N *= 2)
    run();
}

event init (i = 0) {

  /**
  This is the initial Taylor--Green solution for velocity and
  pressure. */
  
  foreach() {
    u.x[] = - cos(2.*pi*x)*sin(2.*pi*y);
    u.y[] =   sin(2.*pi*x)*cos(2.*pi*y);
    p[]   = - (cos(4.*pi*x) + cos(4.*pi*y))/4.;
  }

  /**
  We also need to define the initial centered pressure gradient (this
  improves the accuracy of the initial conditions). */

  foreach()
    foreach_dimension()
      g.x[] = - (p[1] - p[-1])/(2.*Delta);
}

event logfile (i++) {

  /**
  We log the evolution of the maximum divergence and of the total
  kinetic energy. */
  
  scalar div[], ke[];
  foreach() {
    div[] = (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
    ke[] = sq(u.x[]) + sq(u.y[]);
  }
  printf ("%d %d %g %g %g\n", N, i, t, normf(div).max, statsf(ke).sum);
}

event error (t = 2) {

  /**
  At $t=2$ we compute the error on the norm of the velocity. */
  
  scalar e[];
  foreach() {
    double u0 = - cos(2.*pi*x)*sin(2.*pi*y);
    double v0 =   sin(2.*pi*x)*cos(2.*pi*y);
    e[] = norm(u) - sqrt(sq(u0) + sq(v0));
  }
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

/**
For this particular case, the Bell--Collela--Glaz advection scheme
converges at third-order.

~~~gnuplot Accuracy of the solution as a function of the level of refinement
set xlabel 'Spatial resolution'
set ylabel 'Error norms'
set cbrange [1:1]
set logscale
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
     'log' u 1:4 t '|e|_{max}' lc 0, \
     '../taylor-green-all-mach/log' u 1:3 t '|e|_2 (all Mach)', \
     '' u 1:4 t '|e|_{max} (all Mach)'
~~~

The divergence of the centered velocity field is well-behaved.

~~~gnuplot Evolution of the maximum divergence of the centered velocity field
reset
set xlabel 'Time'
set ylabel 'Maximum divergence'
set cbrange [1:1]
set logscale y
set xrange [0:2]
set yrange [1e-4:]
plot '< grep "^32 " out' u 3:4 w l t '32^2', \
     '< grep "^64 " out' u 3:4 w l t '64^2', \
     '< grep "^128 " out' u 3:4 w l t '128^2', \
     '< grep "^256 " out' u 3:4 w l t '256^2', \
     '< grep "^32 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '32^2 (all Mach)',		  \
     '< grep "^64 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '64^2 (all Mach)',		   \
     '< grep "^128 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '128^2 (all Mach)',		   \
     '< grep "^256 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '256^2 (all Mach)'
~~~

By fitting the decrease of the kinetic energy, we get an estimate of
the numerical viscosity.

~~~gnuplot Equivalent Reynolds number as a function of resolution
reset
set xlabel 'Resolution'
set ylabel 'Equivalent Reynolds number'
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

set print "Re-all-mach"
fit [1:] f(x) '< grep "^32 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 32,Re(b)
fit [1:] f(x) '< grep "^64 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 64,Re(b)
fit [1:] f(x) '< grep "^128 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 128,Re(b)
fit [1:] f(x) '< grep "^256 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 256,Re(b)

plot 'Re' w lp t 'centered', 'Re-all-mach' w lp t 'all Mach'
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/reynolds.html)
*/
