/**
# Test dump alternatives
*/
#include "grid/octree.h"
#include "utils.h"
#include "dump2.h"
#include "view.h"

scalar m[];

int main() {
  init_grid (32);
  refine (sq(x - 0.6) + sq(y - 0.4) + sq(z - 0.2) < 0.15 && level < 8);
  foreach()
    m[] = 2.*(x + y + z) + pid();
  dump ("dump_file");
  dump2 ("dump_file2");
  if (pid() == 0)
    system ("cat dump_file2-0* > dump_file2");
  //MPI_Barrier (MPI_COMM_WORLD);
  dump3 ("dump_file3");
  MPI_dump ("dump_file4");
  
  view (fov = 34.2325, quat = {-0.266527,0.362671,0.049675,0.891608},
	tx = -0.000114766, ty = 0.0126711);
  
  restore ("dump_file");
  isosurface ("m", 5);
  box();
  save ("reference.png");
  
  restore ("dump_file2");
  clear ();
  isosurface ("m", 5);
  box();
  save ("dump2.png");
  
  restore ("dump_file3");
  clear ();
  isosurface ("m", 5);
  box();
  save ("dump3.png");

  restore ("dump_file4");
  clear ();
  isosurface ("m", 5);
  box();
  save ("dump4.png");
}
/**
## Did they work? 

This page is run with 6 cores.

![Reference](testd2/reference.png)


![dump2](testd2/dump2.png)


![dump3](testd2/dump3.png)


![dump4](testd2/dump4.png)

*/
