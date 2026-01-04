/**
# Speed of elementary operations on different grids

This code is compiled using either the default tree implementation
or the 'multigrid' regular Cartesian grid implementation.

We use a square, regular, Cartesian grid and vary the resolution from
16^2^ to 2048^2^ (to check for the influence of memory caching). */

#include "utils.h"

scalar a[], b[];

int main (int argc, char * argv[])
{
  size (1.[0]); // dimensionless
  int start, end;
  if (argc == 1) {
    start = 4; end = 11;
  }
  else
    start = end = atoi(argv[1]);
  for (int l = start; l <= end; l++) {
    init_grid (1 << l);
    int nloops, i;

    /**
    We fill `a` with a simple function. */

    foreach()
      a[] = cos(2.*pi*x)*sin(2.*pi*y);

    /**
    We set a number of loops proportional to the number of grid
    points so that times are comparable on all grids. */

    nloops = i = (1 << 25) >> 2*l;

    /**
    Here we compute
    $$
    b = \nabla^2 a
    $$
    using a 5-points Laplacian operator. */
    
    timer start = timer_start();
    while (i--)
      foreach()
	b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[])/sq(Delta);
    printf ("lap %d %g\n", l, 1e9*timer_elapsed (start)/(nloops*(1 << 2*l)));
    
    /**
    Something simpler: the sum of `a` over the entire mesh. */

    scalar s[];
    foreach()
      s[] = a[];
    restriction ({b});
      
    nloops = i = (1 << 25) >> 2*l;
    start = timer_start();
    double sum;
    while (i--) {
      sum = 0.;
      foreach (reduction(+:sum))
	sum += a[];
    }
    printf ("sum %d %g %g\n", l, 1e9*timer_elapsed (start)/(nloops*(1 << 2*l)), sum);

    /**
    And finally the restriction operator. */

    nloops = i = (1 << 25) >> 2*l;
    start = timer_start();
    while (i--)
      restriction ({b});
    printf ("res %d %g %g\n", l, 1e9*timer_elapsed (start)/(nloops*(1 << 2*l)), sum);
  }
}

/**
## Results

This graph shows the speed of the tree implementation for each
operation relative to the speed on the Cartesian mesh. As expected,
the overhead is relatively larger for the simpler operations (e.g. sum
of all elements). This graph is quite sensitive to the exact machine
architecture (cache hierarchy etc...).

~~~gnuplot Relative speed of simple operations on a tree mesh
set xlabel 'Level'
set ylabel 'Cartesian speed / Quadtree speed'
set key top right
set logscale y
plot '< paste out cout | grep lap' u 2:($3/$6) w lp t '5-points Laplacian', \
     '< paste out cout | grep sum' u 2:($3/$7) w lp t 'Sum', \
     '< paste out cout | grep res' u 2:($3/$7) w lp t 'Restriction'         
~~~

The absolute speed for the Laplacian operator on both grid
implementations is shown below. Note that Cartesian meshes are fast!
(i.e. hundreds of million of grid points per second).

~~~gnuplot Absolute speed of the 5-points Laplacian on Cartesian and tree meshes
set ylabel 'nanoseconds per grid point'
set yrange [:]
plot '< grep lap out' u 2:3 w lp t 'Laplacian (Quadtree)', \
     '< grep lap cout' u 2:3 w lp t 'Laplacian (Cartesian)', \
     '< grep sum out' u 2:3 w lp t 'Sum (Quadtree)',   \
     '< grep sum cout' u 2:3 w lp t 'Sum (Cartesian)', \
     '< grep res out' u 2:3 w lp t 'Restriction (Quadtree)',   \
     '< grep res cout' u 2:3 w lp t 'Restriction (Cartesian)'
~~~
*/
