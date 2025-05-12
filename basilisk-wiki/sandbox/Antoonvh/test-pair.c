/**
# A test for the pairwise solver

The solution `s` is split over two domains
   
   ![s1](test-pair/s1.png)
   ![s2](test-pair/s2.png)
 */
#include "poisson-pair.h"
#include "run.h"
scalar s1[] ,s2[];

s1[left] = dirichlet (x);
s1[right] = dirichlet (x);

int main() {
  init_grid (128);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sqrt(sq(x - 0.423) + sq(y - 0.313)) - 0.20789);
  fractions (phi, cs, fs);
  restriction ({cs, fs}); //no auto?
  scalar b1[], b2[];
  foreach() {
    b1[] = b2[] = 0.0;
  }
  face vector k1[], k2[];
  foreach_face() {
    k1.x[] = fs.x[]*kappav1;
    k2.x[] = (1 - fs.x[])*kappav2;
  }
  
  poisson_pair ({s1, s2}, {b1, b2}, k1, k2);
  
  output_ppm (s1, file = "s1.png", n   = 300, min = 0, max = 1);
  output_ppm (s2, file = "s2.png", n   = 300, min = 0, max = 1);
}
