#include "refine_unbalanced.h"

#define BGHOSTS 2

int main (int argc, char * argv[])
{
  X0 = Y0 = -0.5;
  init_grid (2);

  int depth = argc > 1 ? atoi(argv[1]) : 6;
  refine_unbalanced ((level <= depth && x <= -0.25 && y < 0 && y >= -0.25) ||
		     (level <= depth - 1 && y < 0), NULL);
  output_cells (stdout);
  
  scalar s[];
  foreach()
    s[] = 1.;

  // check boundary conditions on leaves
  foreach()
    foreach_neighbor()
      assert (s[] == 1.);
    
  // rebalancing
  int nf = 0;
  foreach(reduction(+:nf))
    nf++;
  int npe;
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  nf = max(1, nf/npe);
  scalar index[];
  z_indexing (index, true);
  foreach()
    fprintf (qerr, "%g %g %g %d\n", x, y, z, min(npe - 1, (int)(index[]/nf)));
}
