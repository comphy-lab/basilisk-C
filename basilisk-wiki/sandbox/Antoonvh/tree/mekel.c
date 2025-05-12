/**
# A Tree in Mekel Park

A poor reconstruction of a tree in TU Delft's Mekel park is shown
[here](realtree.c). Now we take a different approach: The same
point-cloud data is fed to TU Delft's own
[`adtree`](https://github.com/tudelft3d/adtree/blob/master/README.md)
code, and [Meshlab](http://www.meshlab.net/) is then used to convert
the reconstructed tree to a `.stl` file.

![Tree with local curvature](mekel/mekel.mp4)

Much better!
*/
#include "grid/octree.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "curvature.h"
#include "view.h"
  
int main () {
  coord * p = input_stl (fopen ("mekel.stl", "r"));
  coord min, max;
  bounding_box (p, &min, &max);
  //Reposition and rotate
  int j = 0;
  while (p[j].x != nodata) {
    p[j].x -= (min.x + max.x)/2;
    p[j].y -= (min.y + max.y)/2;
    p[j].z -= min.z;
    j++;
  }
  j = 0;
  while (p[j].x != nodata) {
    coord pt = {p[j].x, p[j].y, p[j].z};
    foreach_dimension()
      p[j].x = pt.y;
    j++;
  }
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  init_grid (16);
  size (1.2*maxl);
  X0 = (max.x + min.x)/2. - L0/2;
  Z0 = (max.z + min.z)/2. - L0/2;
  Y0 = (max.y + min.y)/2. - L0/2;
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){3e-4*L0}, 10).nf);
  boundary ({d});
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  scalar f[], kappa[];
  f.prolongation = fraction_refine;
  face vector s[];
  fractions (phi, f, s);
  boundary ({f});
  curvature (f, kappa);
  view (width = 600, height = 900, fov = 17, ty = -0.4);
  for (double th = 0; th < 3*pi; th += 0.05) {
    view (phi = 0.3, theta = th);
    draw_vof ("f", "s", color = "kappa");
    save ("mekel.mp4");
    printf ("Saving frames. Progress: %g %% \n", 100*th/(3*pi)); 
  }
}
