/**
## Packing of the droplets using [packing geometry
function](packing.h) 

![Packing of the
droplets](test_packing2/packing.png)
 */

#include "grid/multigrid.h"
#include "run.h"
#include "fractions.h"
/** Define D and LEVEL */
#define D 1. 
#define LEVEL 9
#include "packing.h"
#include "../../Antoonvh/profile5c.h"
#include "view.h"

scalar f[];

int main()
{
  /** 
## User interface
 */
  i0 = 1.25*D; // Initial distance between bottom to the center of botom
	      //most layer of droplets
  num = 10; // Number of droplets in x-direction
  numy = 15; // Number of droplet layers in y-direction
  iy0 = 1.25*D; // Incremental y-distance
  ix0 = 2.*D; // Incremental x-distance
  pk = 2; 
  double length = 2.*num;
  ar[2] = 0.5 + cell_gap(2, length); //Initialize few cells away from
				     //boundaries
  ar[3] = 1.5 - cell_gap(2, length);
  
  size(length);
  init_grid (1 << LEVEL);
  
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
   set grid
   set xlabel 'f'
   set ylabel 'y/LENGTH'
   set xrange[0:1]
   set yrange[0:1]
   plot 'out' u 2:($1/20) w lp pt 6 ps 0.7 notitle
   ~~~

*/