/**
# Test case for the output of two capsules in Paraview across periodic
boundaries

*/
#define L0 1.
#define RADIUS .125
#define LEVEL 1
#define LAG_LEVEL 4
#define NCAPS 2
#define PARAVIEW_CAPSULES 1

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  TOLERANCE = HUGE;
  foreach_dimension() periodic(left);
  run();
}

event init (i = 0) {
  activate_biconcave_capsule(&CAPS(0), radius = RADIUS, level = LAG_LEVEL, shift={L0/2., 3*L0/7., 3*L0/4.});
  activate_biconcave_capsule(&CAPS(1), radius = RADIUS, level = LAG_LEVEL);

  pv_output_ascii(fp=stderr);
}