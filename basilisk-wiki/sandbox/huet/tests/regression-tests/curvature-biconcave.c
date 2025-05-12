/**
# Membrane curvature test

In this test, we compute the membrane mean and Gaussian curvatures, the 
Laplace-Beltrami operator of the curvature, as well as the bending force
density in the case of a biconcave capsule.

This file tests the functions in [bending-ft.h](../../src/lagrangian_caps/bending-ft.h) and [curvature-ft.h](../../src/lagrangian_caps/curvature-ft.h).

![Shape of a biconcave capsule](curvature-biconcave/rbc.png){ width=30% }
*/
#define L0 4.
#define RADIUS 1
#define LEVEL 1
#define LAG_LEVEL 4
#define REF_CURV 0
#define E_B 1.

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  TOLERANCE = HUGE;
  run();
}

event init (i = 0) {
  activate_biconcave_capsule(&CAPS(0), radius = RADIUS, level = LAG_LEVEL);
}

event output (i = 1) {
  double total_node_area = 0.;
  for(int i=0; i<CAPS(0).nln; i++)
    total_node_area += compute_node_area(&CAPS(0), i);
  double total_triangle_area = 0.;
  for(int i=0; i<CAPS(0).nlt; i++)
    total_triangle_area += CAPS(0).triangles[i].area;

  lagMesh* mesh = &(CAPS(0));
  for(int i=0; i<mesh->nln; i++) {
    double bending_force = cnorm(mesh->nodes[i].lagForce);
    double lbcurv = laplace_beltrami(mesh, i, true);
    double norm = cnorm(mesh->nodes[i].pos);
    double acosarg = (fabs(mesh->nodes[i].pos.y/norm) > 1) ?
      sign(mesh->nodes[i].pos.y)*1. : mesh->nodes[i].pos.y/norm;
    double theta = acos(acosarg);
    double node_area = compute_node_area(&CAPS(0),i);
    fprintf(stderr, "%.4g %.4g %.4g %.4g %.4g %.4g %d\n", theta, mesh->nodes[i].curv,
      mesh->nodes[i].gcurv, lbcurv, bending_force/node_area, node_area, i);
  }

  view(fov = 18, bg = {1,1,1}, theta = pi/6, psi = 0., phi = pi/8, width=1200, height=1200);
  // view(fov = 18, bg = {1,1,1});
  clear();
  draw_lag(mesh = &(CAPS(0)), lw = .5, facets = true);
  save("rbc.png");
}

/**
# Results
## Mean and Gaussian curvatures

~~~gnuplot
set xlabel "Polar angle {/Symbol q}"
set ylabel "Mean and Gaussian curvatures"
set grid
set xtics ('0' 0, 'π/10' pi/10, 0, 'π/5' pi/5, '3π/10' 3*pi/10, '2π/5' 2*pi/5, \
'π/2' pi/2)
set xrange [0:pi/2]

plot 'log' using (($1 < pi/2) ? $1 : pi - $1):2 w p lc "blue" title "Mean curvature", \
'log' using (($1 < pi/2) ? $1 : pi - $1):3 w p lc "orange" title "Gaussian curvature", \
'rbc_exact_curv.csv' using 1:2 w l lc -1 lw 2 dt 2 title "Exact curvatures", \
'rbc_exact_gcurv.csv' using 1:2 w l lc -1 lw 2 dt 2 title ""
~~~

## laplace-Beltrami operator of the mean curvature

~~~gnuplot
set grid
set xtics ('0' 0, 'π/10' pi/10, 0, 'π/5' pi/5, '3π/10' 3*pi/10, '2π/5' 2*pi/5, \
'π/2' pi/2)
set key bottom left reverse Left
set xrange [0:pi/2]

plot 'log' using (($1 < pi/2) ? $1 : pi - $1):($7 > 11 ? $4 : 1/0) w p lc "blue" title "Surface Laplacian of curvature", \
'rbc_exact_lbcurv.csv' using 1:2 w l lc -1 lw 2 dt 2 title "Exact solution"

~~~

## Total nodal bending force density

~~~gnuplot
set xlabel "Polar angle {/Symbol q}"
set ylabel "Bending force density"
set grid
set xtics ('0' 0, 'π/10' pi/10, 0, 'π/5' pi/5, '3π/10' 3*pi/10, '2π/5' 2*pi/5, \
'π/2' pi/2)
set key bottom left reverse Left

plot 'log' using (($1 < pi/2) ? $1 : pi - $1):5 w p lc "blue" title "Bending force density", \
'rbc_exact_bending.csv' using 1:2 w l lc -1 lw 2 dt 2 title "Exact solution"
~~~

*/
