/**
# Distance attenuation

![In 3D far away particles appear smaller](distanceat/nja.png)

*/
#include "grid/multigrid3D.h"
#include "view.h"
vector u[];
#define BVIEW 1
#include "particles.h" // For init_part. function

int main() {
  X0 = Y0 = Z0 = -L0/2.;
  init_grid (4); //a 4 x 4 x 4 grid;
  init_particles_in_cells();
  scatter (loc, s = 45);
  save ("nja.png");
  free (loc);
}