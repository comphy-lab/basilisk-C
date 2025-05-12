/**
# Breaking out of a foreach loop

When using a Cartesian or multigrid, the `break` statement does not directly exit the foreach loop. But is does on a tree-grid, making the behaviour grid dependent. 
*/

#include "grid/multigrid.h"

int main() {
  init_grid (4); // a 4 x 4 grid;
  int j = 0;
  foreach() {
    printf ("This is printed %d time(s)\n", ++j);
    break;
  }
}

/**

## output

The line is printed `N` times as it `break`s only the row iterator loop. See [here](break_out/out).
*/