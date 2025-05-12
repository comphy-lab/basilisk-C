/**
## Hexagonal close packing of droplets using [dense packing geometry
function](close-packing.h) 

![Hexagonal close packing of the
droplets](test_hexapacking/packing.png)
 */

#include "grid/multigrid.h"
#include "run.h"
#include "fractions.h"
/** User interface */
#define D 1. //Diameter of the each droplet
#define LEVEL 9
#define NUM 10 // Number of droplets in x direction 
#define NUMY 11 // Number of droplets in y direction
#define PK 2 // Hexagonal packing is 2 and Square packing is 1
#define I0 0.75*D // Distance between bottom to the center of bottom
		// most droplets
#define GAP 2 //more than 2 to avoid coalescence while using
	      //no-coalescence
#include "close-packing.h"
#include "../../Antoonvh/profile5c.h"
#include "view.h"

scalar f[];

int main()
{
  size(LENGTH);
  init_grid (1 << LEVEL);
  
  fraction (f, packing_geometry(x,y));
  
  dump("dump");
  
  view(tx = -0.499269, ty = -0.403256);
  box();
  draw_vof("f", lc = {0,0.5,1}, lw = 2);
  save("packing.png");

  profile({f});

  run();
  
}

/**
## Results

   ~~~gnuplot Density profile or volume fraction profile 
   set grid
   set xlabel 'f'
   set ylabel 'y/LENGTH'
   plot 'out' u 2:($1/10.4065) w lp pt 6 ps 0.7 notitle
   ~~~
*/
