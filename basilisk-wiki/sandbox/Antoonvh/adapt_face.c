/**
# Adapt_wavelet cant handle face fields

The act of calling "adapt_wavelet()" alters some (ghost) cell values.

~~~gnuplot Locations of faces with "new" values
set size square
plot 'log' w lines, 'out'
~~~
 */

#include "utils.h"
#include "poisson.h"

face vector uf[];

scalar p[];

astats adapt (void) {
  return  adapt_wavelet ({uf.x, uf.y}, (double[]){0.001, 0.001}, 8);
}

double max_div (void) {
   scalar div[];
  foreach() {
    div[] = 0;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] = fabs(div[])/ Delta;
  }
  return statsf(div).max;
}

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2;
  init_grid (N);
  foreach_face(x) 
    uf.x[] = exp(-sq(x) - sq(y));
  TOLERANCE = 1e-9;
  uf.x.refine = refine_face_solenoidal;
  do {
    project (uf, p);
  } while (adapt().nc);
  project (uf, p); //Because "adapt_wavelet" ruins it (see below)
   
  printf ("# Cells before adaptation %ld,\n", grid->tn);
  printf ("# with max divergence: %g\n", max_div());
  
  face vector uf2[];
  foreach_face()
    uf2.x[] = uf.x[];
  
  astats st = adapt();
    
  printf ("# Cells after adaptation: %ld, refined / coarsened: %d / %d\n", grid->tn, st.nf, st.nc);
  printf ("# The max divergence is now %g\n", max_div());
  
  foreach_face() {
    if (uf.x[] != uf2.x[])
      printf ("%g %g\n", x, y);
  }
  output_cells(fp = stderr);
}
