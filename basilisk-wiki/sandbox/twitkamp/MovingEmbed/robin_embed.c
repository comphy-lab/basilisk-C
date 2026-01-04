/**
# Robin boundary condition on embedded boundaries

This test case is identical to basilisk.fr/sandbox/twitkamp/robin_couette.c but the slip boundary condition is imposed on an embedded boundary

# Numerical implementation 

The dirichlet_gradient_x function is tweaked to impose the no-slip boundary condition (dirichlet(.)) "lambda away" from the embed cut-cell.

In fact this just works for any dirichlet(0.) which means it automatically implements the robin boundary condition
*/


#include "grid/multigrid.h"
#include "embed_navier.h"
#include "navier-stokes/centered.h"

scalar un[];
u.x[embed] = (y > 0 ? dirichlet(1.) : dirichlet(0.5));

//no penetration
u.y[embed] = (y > 0 ? dirichlet(0.) : dirichlet(0.));

face vector muv[];
scalar lambdax[], lambday[];
double slip;
int main(){
  
  L0 = 1;
  // center of the channel at y = 0
  origin (-L0/2., -L0/2.);

  mu = muv;

  stokes = true;
  TOLERANCE = 1e-4;
  periodic(left);
  
  N = 64;


  for (double s = 0.; s <= 1.; s += 0.4){
    slip = s;
    run();
  }
}


event init (t = 0) {
  double EPS = L0/N;
  vertex scalar phi[];
  foreach_vertex(){
    phi[] = intersection((L0/4. + EPS/4.) - y, (L0/4. + EPS/4.) + y);
  }

  // I think most of this is not necessary anymore 
  boundary ({phi});
  fractions (phi, cs, fs);
  boundary({cs, fs});
  restriction({cs, fs});
  fractions_cleanup(cs, fs, smin=1e-4);
  boundary({cs, fs});
  restriction({cs, fs});

  /** 
  ## The slip length is an scalar field attribute defined in embed_navier.h
  */

  foreach(){
    lambdax[] = (y > 0. ? 0. : slip);
    lambday[] = 0.;
  }

  u.x.lambda = lambdax;
  u.y.lambda = lambday;

}

event properties (i++)
{
  foreach_face()
    muv.x[] = 1.*fs.x[];
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  fprintf(stderr, "change %g\n", du);
  // #if DOUBLE 
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
  // #endif
}

/** 
Store the velocity as a function of wall position
*/


event profile (t = end)
{
  char file_slip[80];
  sprintf(file_slip, "data_robin_embed_%g.dat", slip);

  FILE *fp_slip = fopen(file_slip, "a");
  foreach(){
    if (cs[] > 0){
      fprintf(fp_slip, "%g %g\n", y, u.x[]);
    }
  }
}
/**
Excellent agreement
## Results
~~~gnuplot Slip velocity as a function of Slip Length
left_bound = -(1.0 + 1.0/64.0)/4.0
right_bound = (1.0 + 1.0/64.0)/4.0

# Add shaded regions
set object 1 rect from graph 0, graph 0 to left_bound, graph 1 fc rgb "gray" fs solid 0.2 behind
set object 2 rect from right_bound, graph 0 to graph 1, graph 1 fc rgb "gray" fs solid 0.2 behind

set key top left
set xlabel "y"
set ylabel "ux"
set grid

plot 'data_robin_embed_0.dat' u 1:2 ps 1 pt 4 lc rgb "red" title 'Sim. data slip length = 0', \
     'data_robin_embed_0.4.dat' u 1:2 ps 1 pt 4 lc rgb "green" title 'Sim. data slip length = 0.4', \
     'data_robin_embed_0.8.dat' u 1:2 ps 1 pt 4 lc rgb "purple" title 'Sim. data slip length = 0.8', \
     '../couette_robin.dat' u 1:2 w l lt 2 lc rgb "black" title 'Expected couette with robin boundary condition'
~~~
*/