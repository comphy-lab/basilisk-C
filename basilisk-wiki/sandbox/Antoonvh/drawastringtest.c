/**
# Obtain expected bahaviour of `draw_string()`?
In Bview, string colours can be a bit dull when combined with an isosurface. To fix this issue, special care is required for the moment.
*/

#include "grid/octree.h"
#include "view.h"
scalar height[];
int main(){
  init_grid(32);
  X0=Z0=Y0=-L0/2.;
  foreach()
    height[]=y;
  boundary({height});
  char title[100];
  // Loop over a range of Height values
  for (double isoval = -0.75; isoval<=0.75; isoval+=0.01){
    sprintf(title,"Isosurface value = %.2g",isoval);
    clear();
    view(fov=30,phi=0.25,theta=0.5,width=640,height=480);
    box(NULL);
    isosurface("height",isoval);
    // Draw a white, magenta and default string.
    glNormal3d (0, 0, 1); // this fixes the problem but should go in `draw_string()`
    draw_string(title,1,40,{1.,1.,1.},3);
    draw_string("Magenta",2,45,{0.99,0.5,0.99});
    draw_string("Default string");
    save("stringtest.mp4");
  }
}
/**
You may view the resulting movie:

![The string colours do not change when the isosurface is drawn](drawastringtest/stringtest.mp4)
*/