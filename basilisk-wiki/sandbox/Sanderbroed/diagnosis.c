/**
# Complex diagnosis 

This page examplifies how to implement complex diagnosis targeted at
vortex-ring simulations. We want to do two things, 

1. Output horizontal slices
2. Analyze temperature slice for the heating evolution

The first is simple as there exist an `sliceXZ()` function if you know
where to look.
 */

#include "grid/octree.h"
#include "../Antoonvh/slicer.h" 
#include "run.h"
#include "view.h"

scalar s[]; 

/**
Instead of computing a solution for `s` as a solution to a
differential equation, its evolution is simply `define`d.
 */

#define FUNC (exp(-sq(x) - sq(y - sin(1.5*t)) - sq(z))*(x - 0.3*sin(t)))

int main() {
   L0 = 5;
   X0 = Y0 = Z0 = -L0/2;
   N = 32;
   run();
}

event init (t = 0) {
  DT = 0.05;
  foreach()
    s[] = FUNC;
  boundary ({s});
}

event advance (i++) {
  dt = dtnext(DT);
  foreach()
    s[] = FUNC;
  boundary ({s});
}

/**
## A bonus movie

![It reminds me of a pumping heart](diagnosis/mov.mp4)
 */

event movie (t += 0.1) {
  box();
  isosurface ("s", -0.1, color = "s", linear = true);
  isosurface ("s",  0.1,  color = "s", linear = true);
  save ("mov.mp4");
}
/**
## Slices

We can take $x-z$ slices and output them for post processing. We
choose the crossection at $yp = -0.4$ and 128 ($= 2^7$) pixels
resolution.

~~~gnuplot The crossection at $t = 3$\
set xlabel 'z axis [px]'
set ylabel 'x axis [px]'
plot 'slice3' matrix with image title 's-field values'
~~~
 */

event slices (t += 1) {
  double yp = -0.4;
  char fname[99];
  sprintf (fname, "slice%g", t);
  sliceXZ (fname, s, yp, 7);
}

/**
## Slice analysis

Rather than post-processing the slices, we may analyze the `s`
field slices with Basilisk C:

We output the surface area of the slice where the `s`-field values are
higher than 0.1, and 0.2. 

~~~gnuplot Areas 
set xlabel 't'
set ylabel  'Area'
set size square
set grid
plot 'areas' u 1:2 w l lw 2 t '>0.1',\
 '' u 1:3 w l lw 2 t '>0.2'
~~~
*/

#define NR_VALS (2) //There are two values to consider.

event areas (t += 0.1) {
  double yp = -0.4;
  static FILE * fp = fopen ("areas", "w");
  double values[NR_VALS] = {0.1, 0.2};     // Values in an array

  //The header:
  if (t == 0) {
    fprintf (fp, "#t");
    for (int i = 0; i < NR_VALS; i++) 
      fprintf (fp, "\t>%g", values[i]);
    fprintf (fp, "\n");
  }
  
  //The values;
  fprintf (fp, "%g", t);
  for (int i = 0; i < NR_VALS; i++) {
    double A = 0; 
    foreach(reduction (+:A)) {
      if (fabs(y - yp) < Delta/2.)                 // Slice is in cell...
	if (interpolate (s, x, yp, z) > values[i]) // Sufficient value
	  A += sq(Delta);                          // Area of slice in sufficient cell
    }
    fprintf (fp, "\t%g", A);
  }
  fprintf (fp, "\n");
}

event stop (t = 10);

