#include "grid/multigrid.h"
#include "run.h"

int LEVEL = 4;

int main()
  {
    init_grid (32);
    //   origin (-L0/2., -L0/2.);

  L0 = 4;  // h 
  

  dimensions (nx=L0,ny=1);
    run();   
  }

event init (i=0)
{

output_cells (stderr);



}

/**
Originalhe
plot [0:4][0:4] "log" u 1:2 w l t "orginal"


*/
