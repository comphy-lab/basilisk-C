/**
# Solar eclipse

The moon may cast a shade on the Earth. For those who missed it on the
12-th or August 2026, we model this phenomenon Bwatch.

![Very convincing. Noting that a flat-Earth eclipse is also possible](eclipse/eclipse.mp4)
*/
   
#include "grid/multigrid3D.h" // We do not need a grid, only three dimensions
#include "bwatch.h"

int main() {
  FILE * fp = popen("ppm2mp4 eclipse.mp4", "w");

  // The sun is modelled as a parallel light source
  default_lights();
  lights[1].dir = (coord){-1, 0.1, -.2};
 
  watch (fov = 5, O = {5, 1, 10}, nx = 800, ny = 500);
  for (t = -0.3; t < 1.1; t+= 0.005) {
    // The background is a dark universe
    sphere (15, {0,0,0}, mat = {.col= {1,1,1}, .col2 = {50,50,50},
				.n1 = {1,0.1,0}, .dull= true});      
    // Earth is a shiny green sphere, its color could be improved.
    sphere (1, {0,0,0}, mat = {{10, 200, 10}});  
    // Shading-casting moon, with realistic tint.
    sphere (.1, {2*cos(t),.1*cos(t),2*sin(t)}, shading = true,
	    mat = {{246, 241, 213}});  
    
    store (fp);
    plain();
  }
  pclose (fp);
}
