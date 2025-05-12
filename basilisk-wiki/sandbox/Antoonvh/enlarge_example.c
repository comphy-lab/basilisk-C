/**
# Horizontal domain extension

Here we examplify the periodic extension of a solution and grid in the two horizontal directions.

It may be usefull for adaptive domain sizing. 
*/
#include "grid/octree.h"
#include "utils.h"
#include "enlarge_domain.h"
#include "view.h"

scalar m[];
int maxlevel = 6; //Maxlevel for first grid.

void output_png(char * fname){
  clear();
  view(fov = 40, ty = -0.25, theta = -pi/1.4,
       phi = pi/5., width = 700, height = 600);
  cells(alpha = Z0);
  cells(n = {1,0,0}, alpha = X0 );
  box(lw = 5);
  isosurface("m", 2);
  save(fname);
}

int main(){
  periodic(left);
#if dimension > 2
  periodic(front);
#endif
  // Original grid and scalar field
  L0 = 2*M_PI;
  X0 = Z0 = -L0/2;
  init_grid( 1 << 3);
  refine(z/10.+sin(x) + 1. > y && level < maxlevel);
  refine(cos(x) + 4. < y && level < maxlevel-1);
  foreach()
    m[] = 0.5*sin(x) + (y) + 0.5*cos(z);
  boundary({m});
  output_png("original.png");
  // 1st Double domain
  enlarge_hor();
  output_png("double.png");
  // Second doubling
  enlarge_hor();
  output_png("quadruple.png");
}

/**
It seems to function on a single core:

![The original grid and solution](enlarge_example/original.png)

![Horizontally extended solution and grid](enlarge_example/double.png)

![Solution and grid after two horizontal extensions](enlarge_example/quadruple.png)


*/