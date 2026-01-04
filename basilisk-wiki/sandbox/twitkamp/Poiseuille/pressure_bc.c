/**
# Poiseuille flow between two horizontal plates. 

I want to simulate a poiseuille flow between two flat embedded plates. 

In this file, I apply a pressure gradient through the boundary condition. 
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

u.n[left] = neumann(0.);
u.n[right] = neumann(0.);

p[left] = dirichlet((double) N);
p[right] = dirichlet(0.);
double EPS;
int main()
{

  size (1. [0]);
  
  origin (0., -L0/2.);

  stokes = true;
  
  TOLERANCE = 1e-5;
  for (int lvl = 4; lvl <= 8; lvl++){
    N = pow(2, lvl);
    run();
  }
}

scalar un[];

event init (t = 0) {

  /**
  The viscosity is unity.*/
   
  mu = fm;
  /**
  The geometry of the embedded boundary is defined as two plates located at $y = \pm (L0 + 2\Delta)/8$. */
  EPS = L0/N;
  solid (cs, fs, intersection((L0/8. + EPS/4.) - y, (L0/8. + EPS/4.) + y));


  /**
   # Embedded boundary conditions

  No slip boundary condition on the embedded boundaries*/
  
  u.x[embed] = y > 0. ? dirichlet(0.) :  dirichlet(0.);
  u.y[embed] = y > 0. ? dirichlet(0.) :  dirichlet(0.);
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-12)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity and error as a function of resolution*/


#define poiseuille(y, N) -N/2. * (y*y - pow(L0/8. + EPS/4., 2))

event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    if (cs[] > 0.) {
      e[] = u.x[] - poiseuille (y, N);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  foreach(){
    if (cs[] > 0.){
      fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	     N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
    }
  }
  
  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("u.x", spread = -1);
  save ("ux.png");

  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e.png");
  
  char file[80];
  sprintf(file, "data_couette_%d", N);
  FILE * fp = fopen(file, "w");
  
  foreach() {
    if (cs[] > 0.){
      fprintf (fp, "%g %g %g %g %g\n",
          y, u.x[], u.y[], p[], e[]);
    }
  }
}

/**
## Results

~~~gnuplot Velocity profile (N = 128)
L0 = 1.
N = 128.
set xlabel 'y'
set ylabel 'u_x'
plate_loc = (L0/8. + (L0/N/4.))
poiseuille(y, N) = - N/2 * (y*y - plate_loc * plate_loc);
set grid
# set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
# set arrow from 0.5, graph 0 to 0.5, graph 1 nohead
plot [-plate_loc:plate_loc][-0.1:1.1]'data_couette_128' u 1:2 t 'numerics', poiseuille(x,N) t 'theory'
~~~


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