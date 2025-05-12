/**
# Test the colourbar

![A color bar?](tcbar/cb.mp4)(loop)
*/

#include "grid/octree.h"
#include "view.h"
#include "colorbar.h"

int main() {
  X0 = Y0 = Z0 = -L0/2;
  init_grid(N);
  scalar height[], f[];

  foreach() {
    double k = 1 [-1];
    height[] = y*k;
    f[] = y*k - 0.4*sin(2*pi*x*k)*cos(3*pi*z*k);
  }
  double min = -.4, max = .4;


  for (double th = 0 ; th <= 2*pi; th += pi/100) {
    view (fov = 30, tx = 0.15, theta = pi/10, phi = pi/20, samples = 4, width = 600, height = 700);
    isosurface("f", 0, color = "height", map = cool_warm, min = min, max = max);
    colourbar (map = cool_warm, min = min, max = max, pos = {.8*sin(th), .8*cos(th)}, 
	       label = height.name, mid = true, draw_box = true, horizontal = false);
    colourbar (map = jet, min = -1, max = 1, pos = {-.8*sin(th), -.8*cos(th)}, 
	       label = "Label", mid = true, draw_box = true, horizontal = true);
    colourbar (map = gray, min = 0, max = 1, pos = {-.6*cos(th), -.2*sin(th)}, 
	       label = "Gray", horizontal = true);
    colourbar (map = gray, min = 0, max = 1, pos = {-.6*cos(th), -.2*sin(th)}, 
	       label = "Gray", horizontal = true);
    colourbar (map = randomap, size = 10+5*sin(2*th), 
	       min = 0, max = 1, pos = {.2*cos(th), .6*sin(th)}, 
	       label = "Growing", horizontal = true);
    cells();
    box();
    save ("cb.mp4");
  }
}
  
