/**
# Intersecting transparant triangles

A challenge for rasterization, easy for ray casting

![](intersect/intersect.mp4)(loop)
 */
#include "grid/multigrid3D.h"
#include "bwatch.h"

int main() {
  // Two triangles
  coord p1[4], p2[4];
  p1[3].x = nodata;
  p2[3].x = nodata;

  p1[0].x = -0.5;
  p1[0].y = -0.5;
  p1[0].z = 0.5;
  
  p1[2].x = -0.5;
  p1[2].y = 0.5;
  p1[2].z = 0.5;
  
  p1[1].x = 0.5;
  p1[1].y = -0.1;
  p1[1].z = -0.5;
  
  p2[0].x = 0.5;
  p2[0].y = 0.5;
  p2[0].z = 0.5;
  
  p2[2].x = 0.5;
  p2[2].y = -0.5;
  p2[2].z = 0.5;
  
  p2[1].x = -0.5;
  p2[1].y = 0.1;
  p2[1].z = -0.5;
  // Movie:
  FILE * fp = popen ("ppm2mp4 intersect.mp4", "w");
  for (double a = 0; a < 2*pi; a += 0.01*pi) {
    watch (fov = 2, O = {2*sin(a), cos(a), 3} );
    sphere (5, mat = {dull = true});
    triangles (p1, mat = {T = 0.4});
    triangles (p2, mat = {T = 0.4, col = {120, 10, 120}});
    store (fp);
    plain();
  }
  pclose (fp);
}
