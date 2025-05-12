/**
# Leve-set test

This test case tries to recover the level-set function from the VOF field.
The test shape is a circle, defined by interface distance (exact level-set). */

#include "utils.h"
#include "fractions.h"
#include "vof_to_LS.h"

int main()
{
  init_grid (8);
  size (1);
  origin (-0.5, -0.5);

  /**
  We initialise the leve-set field *d* and refine the mesh according
  to the error on this field. */
  
  scalar d[],ls[];
  vector gLS[];

  foreach()
    d[] = (sq(0.1) - sq(x-0.1) - sq(y-0.2));
  
  refine(fabsf(d[]) < 0.01 && level < 8);
    
  adapt_wavelet ({d}, (double[]){1e-2}, 9);


  /**
  We initialise a vertex distance field by interpolating the centered
  distance field, and use this to compute VOF fractions.
  Then, level-set and its gradient is computed by *vof_to_LS* function. */
  
  vertex scalar phi[];
  scalar f[];
  
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    
  fractions (phi, f);
  vof_to_LS(f, ls, gLS);

  output_gfs (stdout);
 }