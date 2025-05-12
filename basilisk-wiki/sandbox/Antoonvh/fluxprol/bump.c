/**
# `saint-venant.h` with flux prolongaiton

This is the mostly the `bump.c` code from the [tutorial](/Tutorial).

Changes include:  

* Periodic boundaries 
* Prolongation would trigger an assertion

![Depth](bump/h.mp4)
![Grid](bump/l.mp4)

We can compare the results against A Cartesian solution and the
default Quadtree implementaiton.

~~~gnuplot 
set xlabel 'Time.'
set ylabel 'Depth'
set grid
plot '../bumpC/Cartesian' u 1:2 w l lw 5 lt rgb 'red' t 'min (Cartesian)',		 \
''    u 1:3 w l lw 5 lt rgb 'red' t 'max (Cartesian)',\
 'MyQT' u 1:2 w l lw 3 lt rgb 'green' t 'min (Flux prol)',		 \
''    u 1:3 w l lw 3 lt rgb 'green' t 'max (Flux prol)',\
'../bumpt/Quadtree' u 1:2 w l lw 2 lt rgb 'blue' t 'min (Default Halo Cells)', \
''    u 1:3 w l lw 2 lt rgb 'blue' t 'max (Default Halo Cells)'
~~~

Not very convincing...
 */
#include "myQT.h"
#include "mySV.h"

#define LEVEL 8

int main() {
  periodic (left);     // This should not matter
  periodic (bottom);   // But it does!
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (t = 0) {
  for (scalar s in all)
    s.prolongation = no_prolongation;
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event graphs (i++) {
  char fname[99];
  sprintf (fname, "%s", GRIDNAME);
  static FILE * fp = fopen (fname, "w");
  stats s = statsf (h);
  fprintf (fp, "%g %g %g\n", t, s.min, s.max);
}

event images (t += 4./300.) {
  output_ppm (h, file = "h.mp4", n = 300);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "l.mp4", n = 300, min = 0, max = LEVEL);
}

event end (t = 4); 

/**
## Adaptation

The error-estimator must be overridden. 
 */
event adapt (i++) {
  h.prolongation = refine_linear;
  adapt_wavelet ({h}, (double[]){4e-3}, maxlevel = LEVEL);
  h.prolongation = no_prolongation;
}
