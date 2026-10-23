#include "grid/quadtree.h"
#include "embed.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"

scalar d[];

int maxlevel = 12;
int minlevel = 5;

int main() {
  coord * p = input_xy (fopen ("STFS_logo.gnu", "r"));

  coord pmin, pmax;
  bounding_box (p, &pmin, &pmax);
  double maxl = max (pmax.x - pmin.x, pmax.y - pmin.y);

  size (4.*maxl);
  origin ((pmax.x + pmin.x)/2. - L0/2.,
          (pmax.y + pmin.y)/2. - L0/2.);
  init_grid (1 << minlevel);

  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1e-4*L0}, maxlevel, minlevel).nf);

  solid (cs, fs, -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.);
  fractions_cleanup (cs, fs);

  dump ("restart");
}