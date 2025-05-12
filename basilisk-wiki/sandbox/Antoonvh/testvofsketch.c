/**
# Test for Spherical interface plotting

We plot two spheres, two smooth vof interfaces, and two times the
facets for a vof field.

![There are atleast two issues](testvofsketch/test.mp4)
 */

#include "grid/multigrid3D.h"
#include "bwatch.h"

scalar dummy[], f1[], f2[], f3[], f4[];

int main() {
  FILE * fp = popen ("ppm2mp4 -r 0.5 test.mp4", "w");
  
  X0 = Y0 = Z0 = -L0/2;
  watch (nx = 800, ny = 600);
  
  for (N = 16; N < 257; N *= 2) {
    init_grid (N);
    fraction (f1, (0.12 - sqrt(sq(x - 0.3) + sq(y - 0.2) + sq(z))));
    fraction (f2, (0.12 - sqrt(sq(x) + sq(y - 0.2) + sq(z))));
    fraction (f3, (0.12 - sqrt(sq(x - 0.3) + sq(y + 0.2) + sq(z))));
    fraction (f4, (0.12 - sqrt(sq(x) + sq(y + 0.2) + sq(z))));
    
    restriction ({f1, f3});

    sketch_vof (f3);
    sketch_vof (f1, mat = {ind = 1.1});
    equiplane (f4, vof = true);
    equiplane (f2, vof = true, mat = {ind = 1.1});
    sphere (0.12, {-0.3, -0.2,0}, );
    sphere (0.12, {-0.3,  0.2,0}, mat = {ind = 1.1});
   
    sphere (6, mat = {dull = true});
  
    store (fp);
    
    plain();
  }
  
  pclose (fp);
}
