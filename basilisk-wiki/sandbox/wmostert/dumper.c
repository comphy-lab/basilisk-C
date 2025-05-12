/**
When compiling, using a Makefile I run:

~~~bash
CC='mpicc -D_MPI=27' make dumper.tst
~~~
*/

#include "grid/multigrid3D.h"
#include "utils.h"

int main()
{

  /**
  Initialize the grid. The error occurs only for refinement levels of
  8 and higher. */

  init_grid (1 << 8);

  /**
  Dump with a debug message. */
  
  dump ("test");
  return 0; // return successfully (not 1 which indicates failure)
}
