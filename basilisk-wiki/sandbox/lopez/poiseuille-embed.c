/**
# Poiseuille flow on a slender pipe
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"

uf.n[bottom] = 0.;

#define A 0.21

int main()
{
  origin (-0.5, 0.);
  periodic (right);
  
  stokes = true;
  TOLERANCE = 1e-7;
  
  for (N = 32; N <= 128; N *= 2)
    run();
}

scalar un[];

event init (t = 0) {

  /**
  The gravity vector is aligned with the pipe and viscosity is
  unity. */

  const face vector g[] = {1., 0.};
  a = g;
  mu = fm;

  /**
  The pipe geometry is defined using the levelset function
  $\phi$. A is the radius of the pipe. */  

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (A - y);
  boundary ({phi});
  fractions (phi, cs, fs);

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  restriction ({cm, fm, cs, fs});
    
  /**
  The boundary condition is zero velocity on the embedded boundaries. */

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

  for (scalar s in {u})
    s.third = true;
  
  /**
  We initialize the reference velocity in the fluid domain. */
  
  foreach() {
    u.x[] = (cs[] > 0. ? - 0.25*(sq(A)-sq(y)) : nodata);
    un[] = u.x[];
  }
  boundary({un,u});
}

/**
We check for a stationary solution. */

event logfile (t += 0.1; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/**
We compute the error and display the solution using bview. Eventually
we could plot velocity profiles. */

event profile (t = end) {
  scalar e[];
  foreach() {
    e[] = u.x[] + 0.25*(sq(A)-sq(y));
    //    fprintf (stdout, "%g %g \n", y, u.x[]);
  }
  norm n = normf (e);
  fprintf (ferr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
  draw_vof ("cs", "fs");
  squares ("u.x", linear = true, spread = -1);
  save ("u.x.png");
  dump("dump");
}

/**
![Velocity field](poiseuille-embed/u.x.png)

The method is almost exact.

~~~gnuplot Error as function of resolution
set xlabel 'Resolution'
set ylabel 'Maximum error'
set logscale
set xtics 4,2,128
set ytics format "% .0e"
plot 'log' u 1:4 w lp t ''
~~~
*/
