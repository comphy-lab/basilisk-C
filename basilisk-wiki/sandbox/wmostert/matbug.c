#include "utils.h"

/**
  Check that output_matrix uses the correct y-coordinate.*/

scalar xfield[], yfield[];
int i;

int main(){
  init_grid (256);
  origin (-1., 0.);
  /** Assign two scalar fields to the values we know x and y to take on the domain. */
  foreach() {
    xfield[] = x;
    yfield[] = y;
  }
  /** And plot these fields on output_matrix. If all is well, then a contour map of each field should show x-axis on (-1, 0) and y-axis on (0, 1), and the scalar field xfield should range vary over (-1, 0) and yfield over (0, 1). */
  
  static FILE * fpx = fopen ("xf.dat", "w");
  static FILE * fpy = fopen ("yf.dat", "w");
  output_matrix (xfield, fpx, 256, false);
  output_matrix (yfield, fpy, 256, false);
  return 0;
}

/** We plot the resulting binary files.
~~~gnuplot Test of output_matrix
set term pngcairo enhanced size 640,1120 font ",8"
set size ratio -1
set output "output_matrix.png"
set multiplot layout 1,2 scale 1,1
set pm3d map
splot 'xf.dat' binary u 2:1:3
splot 'yf.dat' binary u 2:1:3
unset multiplot
~~~

The colour on the left ranges from -1 to 0, and on the right from 0 to 1, as expected. However, on both plots the axes extend from -1, 0. This is incorrect. The expected behaviour is that the y-axis should extend from 0 to 1.

*/
