/**
# Refinement of axisymmetric metric and face fields

Checks that refinement of axisymmetric metric and face fields behaves
properly. */

#include "axi.h"
#include "navier-stokes/centered.h"

double R0 = 0.1;
int LEVEL = 6;

uf.t[top] = dirichlet(0);

int main()
{
  L0 = 2;
  N = 16;
  run();
}

event init (t = 0)
{
  foreach_face()
    uf.x[] = 0.;
  foreach_face (x)
    uf.x[] = fm.x[]*(sq(L0) - y*y);
  
  refine (level < LEVEL && sq(x - 2.) + sq(y) < sq(1.5*R0));

  output_cells (stdout);
  
  foreach()
    assert (fabs (cm[] - y) <= 1e-20);
  
  foreach_face(y) {    
    // fprintf (stderr, "%g %g %g %g\n", x, y, fm.y[], fm.y[] - y);
    assert (fabs (fm.y[] - y) <= 1e-20);
  }
  
  foreach_face(x) {
    // fprintf (stderr, "%g %g %g %g\n", x, y, fm.x[], fm.x[] - y);
    assert (fabs (fm.x[] - y) <= 1e-20);
  }

  foreach_face (x)
    fprintf (stderr, "%g %g %g %g\n", y, uf.x[], cm[], uf.x[]/cm[]);
}

/**
~~~gnuplot Exact and approximate solutions after refinement.
set xlabel 'y'
plot 'log' u 1:($2/$3) t 'uf.x/cm', (4. - x*x) t 'uf.x/cm exact', \
     '' u 1:2 t 'uf.x', x*(4. - x*x) t 'uf.x (exact)', \
     '' u 1:3 t 'cm', x t 'cm (exact)'
~~~
*/
