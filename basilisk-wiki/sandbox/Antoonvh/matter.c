/**
We load an [stl file](https://www.thingiverse.com/thing:247219) from
the [Collectible Mountain
collection](https://www.thingiverse.com/Shapespeare/collections/10km-1100000-collectible-mountains),
by [Sphapespeare/Withney
Potter](https://www.thingiverse.com/Shapespeare/about) (CC 4.0 BY NC).
*/

#include "grid/octree.h"
#include "distance.h"
#include "fractions.h"
#include "domain_dohicky.h"
#include "view.h"

int LEVEL = 9;
void fraction_from_stl (scalar f, FILE * fp, double eps, int maxlevel) {
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double a = 0;
  foreach_dimension()
    if (max.x - min.x > a)
      a = max.x - min.x;
  L0 = a;
  X0 = Y0 = -L0/2;
  Z0 = min.z;
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel, 5).nf);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  fractions (phi, f);
}

scalar f[];

int main () {
  init_grid (64);
  f.refine = f.prolongation = fraction_refine;
  system ("wget https://www.thingiverse.com/thing:247219/zip");
  system ("unzip zip");
  system ("mv Matterhorn_10km_Collectible_Mountain/matterhorn10k.stl matter.stl");
  FILE * fp = fopen ("matter.stl", "r");
  fraction_from_stl (f, fp, 4e-4, LEVEL);
  fclose (fp);
  view (theta = pi/2, tx = 0.35);
  draw_vof("f");
    box();
  save ("f.png");

  /**
![This is not the optimal perspective](matter/f.png)

We wish to swap directions such that the current vertical (`z`) will
become `y`: x -> z -> y -> x
   */
  swap_direction ({f});
  boundary ({f});
  box();
  draw_vof("f");
  save ("f2.png");
  /**
![Well done `swap_direction()'!](matter/f2.png)
   */
}
