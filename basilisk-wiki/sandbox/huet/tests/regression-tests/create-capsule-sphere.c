/**
# Capsule creation test

In this case we simply create a spherical capsule of radius 1 and ensure that
all the membrane nodes are at the right location. We don't test the connectivity
as it should be a necessary condition for the nodes position to be correct.

This file essentially tests most functions in geometry-ft.h and 
common-shapes-ft.h.
*/

#define RADIUS 1.
#define NCAPS 1
#define LAG_LEVEL 4

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"

int main() {
    origin(-.5*L0, -.5*L0, -.5*L0);
    foreach_dimension() periodic(left);
    run();
}

event init (i = 0) {
    activate_spherical_capsule(&CAPS(0), radius=RADIUS, level=LAG_LEVEL, 
        shift={L0/2., 3*L0/7. -2*L0/5});
    for(int i=0; i<CAPS(0).nln; i++) {
        foreach_dimension() fprintf(stderr, "%g ", CAPS(0).nodes[i].pos.x);
        fprintf(stderr, "\n");
    }
    exit(0);
}