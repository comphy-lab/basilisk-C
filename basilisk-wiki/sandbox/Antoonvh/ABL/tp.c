/**
# Generate periodic Topography 

![Perlin's noise](tp/h.mp4)

are they really periodic?

![Yes, there are no issues when the last
 topography is copied four times.](tp/h.png)
 */

#include "grid/octree.h"
#include "perlin.h"
#include "view.h"
#include "../domain_dohicky.h"

scalar cs[];
face vector fs[];
vertex scalar phi[];

int main() {
  periodic (left);
  periodic (back);
  L0 = 11;
  X0= Y0 = Z0 = -L0/2;
  init_grid (32);
  refine (fabs(y) < 1 && level < 7);

  view (theta = 0.3, phi = 0.4);
  for (int g = 0; g < 35; g++) {
    srand (g);
    
    int nx = 5;
    int ny = 3;
    init_perlin (nx, ny);
    foreach_vertex()
      phi[] = y - 2*(perlin (x, z, nx, ny));
    free (gradp);
    nx = 7;
    ny = 7;
    init_perlin (nx, ny);
    foreach_vertex()
      phi[] += 0.75*perlin (x, z, nx, ny);
    free (gradp);
    nx = 17;
    ny = 17;
    init_perlin (nx, ny);
    foreach_vertex()
      phi[] += 0.2*perlin (x, z, nx, ny);
    free (gradp);
    nx = 31;
    ny = 31;
    init_perlin (nx, ny);
    foreach_vertex()
      phi[] += 0.05*perlin (x, z, nx, ny);
    free (gradp);

    scalar h[];
    foreach()
      h[] = y;
    boundary ({phi});
    fractions (phi, cs, fs);
    draw_vof ("cs", "fs", color = "h", min = -1, max = 1);
    save ("h.mp4", opt = "-r 2");
  }
  enlarge_hor({cs}); //Bview could do `translate()`...
  foreach() {
    if (y > Y0 + L0/3)
      cs[] = 1;
  }
  draw_vof ("cs");
  save ("h.png");
}
