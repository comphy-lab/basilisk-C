// generates the dump file read by restore-rectangle.c

#include "grid/multigrid.h"
#include "utils.h"

int main()
{
  int depth = 6;
  dimensions (nx = 2, ny = 1);
  size (1[0]);
  origin (-0.5, -0.5, -0.5);
  init_grid (1 << depth);
  
  scalar s[];
  foreach()
    s[] = sin(x)*cos(y) + 2.*x + y;

  output_cells (stdout);
  dump (file = "restore-rectangle.dump", list = {s});  
}
