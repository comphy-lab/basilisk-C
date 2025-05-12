/**
# Couette flow between rotating cylinders

We test embedded boundaries by solving the (Stokes) Couette flow
between two rotating cylinders. If SLIP is defined, the outer cylinder
is free of viscous stresses. */

#include "grid/multigrid.h"
#include "src/embed_C.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define SLIP 1
int main()
{
  size (1. [0]);
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);
  
  stokes = true;
  TOLERANCE = 1e-5;
  
#if !SLIP  
  int NF = 32;
#else
  int NF = 32;
#endif
  
  for (N = 16; N <= NF; N *= 2)
    run();
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {

  /**
  The viscosity is unity. */
  
  mu = fm;

  /**
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = difference (sq(0.5) - sq(x) - sq(y),
			sq(0.25) - sq(x) - sq(y));
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
 
#if SLIP
  u.n[embed] = (x*x + y*y > 0.14 ? navier (0.) : dirichlet (-y));
  u.t[embed] = (x*x + y*y > 0.14 ? navier (0.) : dirichlet (x));
#else
   u.n[embed] = dirichlet (x*x + y*y > 0.14 ? 0. : -y);
   u.t[embed] = dirichlet (x*x + y*y > 0.14 ? 0. :   x);
#endif

  /**
  We initialize the reference velocity field. */
  
  foreach()
    un[] = u.y[];
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/**
We compute error norms and display the angular velocity, pressure and
error fields using bview. */

#if SLIP
#define powerlaw(r) (sq(0.25)/(sq(0.25) + sq(0.5))*(sq(0.5)/r + r))
#else
#define powerlaw(r) (r*(sq(0.5/r) - 1.)/(sq(0.5/0.25) - 1.))
#endif

event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (cs[] > 0.) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      e[] = utheta[] - powerlaw (r);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  fprintf (ferr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("utheta", spread = -1);
  save ("utheta.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e.png");

  if (N == 32)
    foreach() {
      double theta = atan2(y, x), r = sqrt(x*x + y*y);
      fprintf (stdout, "%g %g %g %g %g %g %g\n",
	       r, theta, u.x[], u.y[], p[], utheta[], e[]);
    }
}

/**
## Results

![Angular velocity](slip/utheta.png)

![Pressure field](slip/p.png)

![Error field](slip/e.png)

~~~gnuplot Velocity profile (N = 32)
set xlabel 'r'
set ylabel 'u_theta'
powerlaw(r)=(0.25**2/(0.25**2 + 0.5**2)*(0.5**2/r + r))
set grid
set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
set arrow from 0.5, graph 0 to 0.5, graph 1 nohead
plot [0.2:0.55][-0.05:0.35] 'out' u 1:6 t 'numerics',powerlaw(x) t 'theory'
~~~


## See also

* [Wannier flow between rotating excentric cylinders](wannier.c)
*/
