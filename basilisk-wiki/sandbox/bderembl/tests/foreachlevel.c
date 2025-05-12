
/**

# foreach_level with and without MPI

As illustrated below, the foreach_level does not seem to behave the same way with
MPI and without MPI.

In the code I use foreach_level(1)

   witout MPI, I get 4 cells

0.25 	 0.25 	 0.25
0.25 	 0.75 	 0.25
0.75 	 0.25 	 0.75
0.75 	 0.75 	 0.75

which is what I expect for a foreach_level(1)

however with 

4 processors, I get

0.625 	 0.625 	 0.625
0.625 	 0.875 	 0.625
0.875 	 0.625 	 0.875
0.875 	 0.875 	 0.875
0.125 	 0.125 	 0.125
0.125 	 0.375 	 0.125
0.375 	 0.125 	 0.375
0.375 	 0.375 	 0.375
0.125 	 0.625 	 0.125
0.125 	 0.875 	 0.125
0.375 	 0.625 	 0.375
0.375 	 0.875 	 0.375
0.625 	 0.125 	 0.625
0.625 	 0.375 	 0.625
0.875 	 0.125 	 0.875
0.875 	 0.375 	 0.875

which seems to correspond to level 2?

Or should it be understood level(1) per processor?
*/

#include "grid/multigrid.h"
#include "run.h"


scalar omega[];

int main() {
  N = 8;
  init_grid (N);

  foreach() {
    omega[] = x;
  }

  restriction({omega});

  foreach_level(1) {
    printf("%g \t %g \t %g\n",x,y,omega[]);
  }

}
