/**
# `bwatch`; a Basilisk-based ray-casting rendering system

![](bwatcht/nja.png)

*/
#include "grid/octree.h"
#include "bwatch.h"
#include "distance.h"
#include "curvature.h"

int main() {
  X0 = Y0 = Z0 = -L0/2 ;
  N = 64;
  init_grid (N);
  scalar f[];
  coord * p = input_stl (fopen ("bwatcht.stl", "r"));
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  size (2.*maxl);
  X0 = (max.x + min.x)/2. - L0/2;
  Z0 = (max.z + min.z)/2. - L0/2;
  Y0 = (max.y + min.y)/2. - L0/2;
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1e-4*L0}, 8).nf);
  vertex scalar phi[];  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	      d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f);

  scalar s[], kappa[];
  curvature (f, kappa);
  foreach()
    s[] = x + y + z;

  watch (nx = 1024, ny = 768, fov = 1.2*L0,
	 poi = (coord){X0 + L0/2., Y0 + L0/2.,Z0 + L0/2.},
	 O = (coord){X0 + L0/2., Y0 + L0/2., Z0 + 5*L0});
  quadriangles (d, -0.3*L0, mat = {.min = -L0/2., .max = 0.1});
   //fix me: The lattice must be with "HQ"
  lattice (-0.29*L0, HQ = true);
 
  sketch_vof (f, mat = {.s = kappa, .map = cool_warm, .min = -0.05, .max = 0.5, .R = 0.05});
  store (fopen ("nja.ppm", "w"));
  plain();
  system ("convert nja.ppm -resize 80% nja.png");
}



