/**
# Capsule creation test

In this case we draw a triangulated capsule on top of an Eulerian field
to test [view-ft.h](../../src/lagrangian_caps/view-ft.h).
*/

#define L0 1.
#define RADIUS .25
#define NCAPS 1
#define LAG_LEVEL 4

#include "grid/octree.h"
// fixme: capsules should be independent from the solver
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

int main() {
  origin (-.5*L0, -.5*L0, -.5*L0);
  run();
}

event init (i = 0) {
  activate_spherical_capsule (&CAPS(0), radius = RADIUS, level = LAG_LEVEL);
  foreach()
    u.x[] = sin(2*pi*x*y/L0);

  view (fov = 18.9, bg = {1,1,1}, width=300, height=300);
  clear();
  draw_lags (lw = .5, ns = 1.5, facets = true);
  squares ("u.x", n = {0,0,1});
  cells (n = {0,0,1});
  save ("sphere.png");
}

/**
## Result

![Generated image](draw-sphere/sphere.png)

![Target image](http://basilisk.fr/sandbox/huet/tests/regression-tests/draw-sphere/sphere-target.png)
*/
