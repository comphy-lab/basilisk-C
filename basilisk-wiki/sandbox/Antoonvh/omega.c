/**
# Dual-grid AMR solver for the $\omega-\psi$ formulation of the Navier-Stokes equations

We split the advection of vorticty ($\omega$) and the Poisson-equation
solving for the stream function ($\psi$) on separated adaptive
grids. This uses the linked object file similar to the [master-slave](/src/master.h) design.

Results:

![The vorticity field ($\omega$) represented on the fine mesh (left)
 and teal streamlines ($\psi = \text{Const.}$) represented on the
 coarse mesh (right) for a case of two merging
 vortices](https://www.antoonvanhooft.nl/media/omgpsi_res.mp4)
 */
#include "master-omgpsi.h"

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2;
  DT = 0.01;
  run();
}

event init (t = 0) {
  foreach()
    omega[] = exp(-sq(x-1) - sq(y)) + exp(-sq(x + 1) - sq(y));
}

#if MOVIE
#include "view.h"
event movs (t += 0.25) {
  view (fov = 20);
  squares ("omega", map = blue_white_red, min = -0.75, max = 0.75, linear = true);
  cells();
  save ("omega.mp4");
  save ("omega.png");
  slave_level();
}
#endif

event adapt (i++, last) {
  adapt_wavelet ({omega, psi}, (double[]){0.01, 0.01}, 8);
}

event stop (t = 50);
