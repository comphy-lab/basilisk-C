/**
# MWE for mirror

![Dark mirror image](mirror_mwe/mirror.png)

*/

#include "fractions.h"
#include "view.h"
#define circle(x, y, R) (sq(y) + sq(x) - R) // circular initial interface

int main() {
  size (3.);
  origin (0., 0.); 
	N = 1 << 6 ;
	init_grid(N);

  scalar f[];
  fraction (f, circle(x, y, 1.));

  /**
  We prepare the image with view(). */

  view (fov = 30, width = 320, height = 320, samples = 1);
  clear();
  squares ("f", min = 0., max = 1.);
  
  /**
  The mirror image is done with the mirror() function. */
  
  mirror({0,-1,0}) {
    squares ("f", min = 0., max = 1.);
  }
  
  /**
  We save a .png image. */
  
  save ("mirror.png");
}
