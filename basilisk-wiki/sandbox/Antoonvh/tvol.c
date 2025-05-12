/**
# Test for volumetric rendering

It does what it should

![Grid and Gaussian smoke](tvol/vol.mp4) 
*/

#include "grid/multigrid3D.h"
#include "bwatch.h"
scalar dummy[];
scalar s[];
int main() {
  X0 = Y0 = Z0 = -L0/2;
  init_grid(N);
  foreach() 
    s[] = exp(-(sq(x/1.2)+sq(y*1.2)+sq(z))*20);
  quadriangles (s, alpha = -0.4);
  volume (s, cols = true, sc = 0.1, mval = 0.1, shading = 1);
  store (fopen("smk.ppm", "w"));
  
  
  //return 0;
  
  FILE * fp = popen ("ppm2mp4 vol.mp4", "w");

  // Adjust transparancy
  for (double a = 0.01; a < 0.1; a *= 1.05) {
    quadriangles (s, alpha = -0.4);
    volume (s, cols = true, shading = 1, sc = a, mval = 0.1);
    store (fp);
    plain();
  }
  // move camera
  for (double a = 0.; a <= 2*pi; a += pi/50) {
    coord Oc = {3*sin(a), sin(2*a), 3*cos(a)};
    watch (O = Oc);
    quadriangles (s, alpha = -0.4);
    volume (s, sc = 0.1, cols = true, shading = 1, mval = 0.1);
    store (fp);
    plain();
  }
  fflush (fp);
  pclose (fp);
}
