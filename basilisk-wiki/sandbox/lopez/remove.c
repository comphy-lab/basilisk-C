/**
# Test for withdrawal of droplets */

//#include "grid/cartesian.h"
#include "run.h"
#include "tag.h"
#include "src/droplet_stat.h"

scalar c[];

int main()
{
  L0 = 0.5 [0];
  origin (-L0, -L0);
  init_grid (64);
  run ();
}

/**
A set of circles of decreasing radius are distributed in the
computational domain.
  
A proper initialization would need the Vofi library to compute the
correct volume fractions but we do a rough initizalization setting
to 1 these circles. */

event init (i = 0) {  
  foreach()
    c[] = 0.;
  
  int nx = 3, ny = 2;
  double radius = 0.1;
  for (int j = 1; j <= ny; j++) 
    for (int i = 1; i <= nx; i++) {
      double xc = -0.5 + (double) i/(nx + 1); 
      double yc = -0.5 + (double) j/(ny + 1);
      radius *= 0.8;
      foreach()
	if(sq(x-xc)+sq(y-yc) < sq(radius))
	  c[] = 1.;
    }
}

event drop_remove (t = 0) {
  droplet_remove (c, 3, true); 
  output_field ({c}, stderr);    
}

/**
~~~gnuplot The three largest droplets remain  (three removed)
set pm3d map 
set size ratio -1
splot 'log' u 1:2:3
~~~
 */
