/**
# Spiral on a sphere

Since [Desmos](www.desmos.com), does not have proper support for spherical coordinates, we (must) "do" it here.

We study the intersection curve between the planes given by,
$$r = 1,$$
$$\phi = \theta$$

where $r$ is the radial coordinate, $\phi$ is the azimuth angle and
$\theta$ is the polar angle. We use bview to visualize.

![A white sphere ($r = 1$) intersects with the teal $\theta = \phi$ surface along a grey curve. The Red, yellow and green arrows indicate the $x,y$ and $z$ directions, respectively](spiral_on_sphere/spiral.png){width=350px}

![A movie with changing perspectives rendered with bview.](spiral_on_sphere/helix.mp4)

 */
#include "grid/multigrid3D.h"
#include "view.h"
#include "fractions.h"
#include "some_primitives.h"
#define cross cross2
#include "bwatch.h"
// Define spherical coordinates away from the origin

#define RAD (sqrt(sq(x) + sq(y) + sq(z)))
#define THETA (acos(z/RAD))
#define PHI (sign(y)*acos(x/sqrt(sq(x) + sq(y))))

coord curve (double theta) {
  return (coord){cos(theta)*sin(theta), sin(theta)*sin(theta), cos(theta)};
}

int main() {
  L0 = 2.5;
  X0 = Y0 = Z0 = -L0/2.0123456789;
  init_grid (128);
  scalar f1[], f2[];

  view(quat = {0.284, 0.587, 0.674, 0.347});
  fraction (f1, RAD - 1);
  fraction (f2, THETA - PHI);
  restriction ({f1, f2});
  draw_vof ("f1");
  draw_vof ("f2", fc = {0.0, 0.7, 0.7});
  draw_axis_cross (len = 1.2);
  draw_paramterization (0.03, 0, pi,  curve, 50, fc = {0.5, 0.5, 0.5});
  save ("spiral.png");
  /**
We also attempt to visualize this using bwatch. 
   */
  // A scalar field to visualize the curve
  scalar curve_s[];
  foreach() {
    double dist = HUGE;
    for (double gamma = 0; gamma < pi; gamma += 0.01) {
      coord a = curve(gamma);
      if (sq(x - a.x) + sq(y - a.y) + sq(z - a.z) < dist)
	dist = sq(x - a.x) + sq(y - a.y) + sq(z - a.z);
    }
    curve_s[] = exp(-dist);
  }
  restriction ({curve_s});
  
  watch (O = {3,3,1}, up = {1e-6, 1e-6, 1}, fov = 4);
  lights[1].dir = (coord){-1, -1, -3};            
  equiplane (f2, 0, vof = true, mat = {.col = {1, 180, 180}, .T = 0.3}); // phi = theta surface
  equiplane (curve_s, 0.99, mat = {.col = {158, 158, 158}});             // The curve
  sphere (R = 1, mat = {.col = {200, 200, 200}, .T = 0.4, .R = 0});              // r = 1 surface
  sphere (R = 20, mat = {.col = {70,130,180}, .dull = true});          // Background sphere to catch the rays that pass transparent object.
  
  FILE * mp = popen ("ppm2mp4 helix.mp4", "w");

  for (double t = 0; t < pi/2.; t += 0.01) {
    printf ("%g\n", t);
    watch (O = {3*cos(t), 3*sin(t), sin(4*t)});
    store (mp);
  }
  pclose (mp);
}
