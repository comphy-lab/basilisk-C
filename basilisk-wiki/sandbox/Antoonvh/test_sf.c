/**
# Test the structure function function

It seems to be off by a factor of two
*/
#include "grid/octree.h"
#include "structure_function.h"

int main() {
  init_grid (128);
  vector v[];
  
  foreach() {
    v.x[] = x;
    v.y[] = 0.;
    v.z[] = 0.;
  }

  long2structure (stdout, v, 5000, L0/1.5, 50);
   /**
  ## $\mathbf{v} = \{x, 0\}$
  
  $\delta v_{//} (l) = l cos(\theta)^2$
  
  average $\delta v^2$ over $\theta$ 
  
  $S_2(l) = \frac{3}{8}l^2$
  
  ~~~gnuplot
  plot 'out' u 1:2 t 'data', 0.375*x**2 t 'analytical'
  ~~~
  */
  foreach() {
    v.x[] = 0.;
    v.y[] = x;
  }
  
  long2structure (stderr, v, 5000, L0/1.5, 50);
  
  
    /**
## $\mathbf{v} = \{0, x\}$
  
  $\delta v_{//} (l) = l sin(\theta)cos(\theta)$
  
  average $\delta v^2$ over $\theta$ 
  
  $S_2(l) = \frac{1}{8}l^2$
  ~~~gnuplot
  plot 'log' u 1:2 t 'data', 0.125*x**2 t 'analytical'
   ~~~
  */
}
