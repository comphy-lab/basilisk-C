/**
# 'Compressing' an image using adapt_wavelet

*/


#include "grid/quadtree.h"
#include "input.h"
#include "view.h"
#include "./icon.h"



scalar g[];

int main(){
  init_grid(256);
  input_icon(g);
  view (tx = -0.5, ty = -0.5, fov = 20, width = 1000, height = 1000);
  cells();
  squares("g", min = 1, max = 0, map = gray);
  save("start.png");
  while (adapt_wavelet ({g}, (double []){0.2}, maxlevel = 8).nc);
  cells();
  squares("g", min = 1, max = 0, map = gray);
  save("end.png");
}

/**
<img src="wavelet_comp/start.png" alt="Before 'compression'" width="400"/>
<img src="wavelet_comp/end.png" alt="After 'compression'" width="400"/>
*/