/**
# Couette flow between two horizontal plates. 

This test is heavily inspired by https://basilisk.fr/src/test/couette.c.
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"


int main()
{

  size (1. [0]);
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);

  periodic(right);
  
  stokes = true;
  
  TOLERANCE = 1e-5;
  for (int lvl = 4; lvl <= 8; lvl++){
    N = pow(2, lvl);
    run();
  }
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {

  /**
  The viscosity is unity. */
  
  mu = fm;

  /**
  The geometry of the embedded boundary is defined as two plates located at $y = \pm (L0 + \Delta)/4$. */
  double EPS = L0/N;
  solid (cs, fs, intersection((L0/4. + EPS/4.) - y, (L0/4. + EPS/4.) + y));


  /**
  # Embedded boundary conditions

  The top plate moves with a tangential velocity $u_x = 1$ and the bottom plate with a tangential velocity $u_x = -1$.
 */

  u.x[embed] = y > 0. ? dirichlet(1.) :  dirichlet(-1.);
  u.y[embed] = y > 0. ? dirichlet(0.) :  dirichlet(0.);
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 100000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-7)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity, pressure and
error fields using bview. */


#define couette(y, N) (y / (L0/4. + L0/N/4.))

event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    if (cs[] > 0.) {
      e[] = u.x[] - couette (y, N);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
  if (N == 64){
    char file[80];
    sprintf(file, "data_couette_%d", N);
    FILE * fp = fopen(file, "w");
    foreach() {
      fprintf (fp, "%g %g %g %g %g\n",
          y, u.x[], u.y[], p[], e[]);
    }
  }
}

/**
## Results

~~~gnuplot Velocity profile (N = 64)
L0 = 1.
N = 64.
set xlabel 'y'
set ylabel 'u_x'
plate_loc = (L0/4. + (L0/N/4.))
couette(y, N)=(y / plate_loc)
set grid
# set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
# set arrow from 0.5, graph 0 to 0.5, graph 1 nohead
plot [-plate_loc:plate_loc][-1:1]'data_couette_64' u 1:2 t 'numerics', couette(x,N) t 'theory'
~~~

Error goes down. Slow convergence rate maybe due to already good results even for low N (simplicity of the setup).

The error is dependent on the $du$ condition set in the logfile event...

~~~gnuplot Error convergence
unset arrow
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x) = a2 + b2*x
fit f2(x) '' u (log($1)):(log($2)) via a2,b2
set xlabel 'Resolution'
set logscale
set xtics 8,2,1024
set ytics format "% .0e"
set grid ytics
set cbrange [1:2]
set xrange [8:512]
set ylabel 'Error'
set yrange [*:*]
set key top right
plot '' u 1:4 pt 6 t 'max', exp(f(log(x))) t ftitle(a,b), \
     '' u 1:2 t 'avg', exp(f2(log(x))) t ftitle(a2,b2)
~~~
*/
