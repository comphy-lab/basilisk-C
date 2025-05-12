/**
## Hexagonal close packing of droplets using [dense packing geometry
function]packing_geometry1.h) 

![Hexagonal close packing of the
droplets](test_geometry/packing.png)
 */

#include "grid/multigrid.h"
#include "run.h"
#include "fractions.h"
/** User interface */
#define D 1. //Diameter of the each droplet
#define LEVEL 8
#define NUM 10 // Number of droplets in x direction 
#define NUMY 10 // Number of droplets in y direction
#define DIST (0.2 + D)
#define PK 2 // Hexagonal packing is 2 and Square packing is 1
#define I0 1.25*D // Distance between bottom to the center of bottom
		// most droplets
#define theta_pac 0. //0 for 0 degrees rotaion, 0.523599 for 30
			   //degrees rotation
#define rand_magn 5 //5 represents max(0.05*D) amplitude random
		    //perturbation, for no perturbation set it to 0
#include "close_packing.h"
#include "../../Antoonvh/profile5c.h"
#include "view.h"

scalar f[];

int main()
{
  size(LENGTH);
  init_grid (1 << LEVEL);
  rval(rval1, rval2, NUM, NUMY, rand_magn);
  fraction (f, packing_geometry(x,y));
  
  dump("dump");
  
  view(fov = 21.2684, tx = -0.492457, ty = -0.458159);
  box();
  draw_vof("f", lc = {1,0,0}, lw = 2);
  save("packing.png");

  profile({f});
  run();
  
}

/**
## Results

   ~~~gnuplot Density profile or volume fraction profile 
   reset
   set grid
   set xlabel 'f'
   set ylabel 'y'
   plot 'out' u 2:1 w lp pt 6 ps 0.7 notitle
   ~~~

*/
