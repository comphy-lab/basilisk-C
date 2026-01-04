#include "utils.h"
#include "output.h"
#include "refine_unbalanced.h"
#include "check_restriction.h"

#define BGHOSTS 2

void image()
{
  // const scalar pid[] = pid(); // fixme: this should work
  scalar pid[];
  foreach()
    pid[] = pid();
  output_ppm (pid, n = 128, min = 0, max = npe() - 1);
}

void set_unused_pids()
{
  scalar pid[];
  foreach_cell()
    pid[] = is_local(cell) ? pid() : -1;
  for (int l = 0; l <= depth(); l++)
    boundary_iterate (restriction, {pid}, l);
  foreach_cell()
    if (pid[] < 0)
      cell.pid = npe();
}

void partition (const char * prog)
{
  // generates reference solution (in ref)
  init_grid (1);
  
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
  
  refine_unbalanced (level < 5, NULL);

  refine_unbalanced (level <= 9 && sq(x - 0.5) + sq(y - 0.5) < sq(0.05), NULL);
  
  mpi_partitioning();
  set_unused_pids();

  char name[80];
  sprintf (name, "ref-%d", pid());
  FILE * fp = fopen (name, "w");
  debug_mpi (fp);
  fclose (fp);

  MPI_Barrier (MPI_COMM_WORLD);
  if (pid() == 0) {
    sprintf (name, "cat ref-* > ref && cp -f ref %s.ref", prog);
    system (name);
  }
}

int main (int argc, char * argv[])
{
  partition (argv[0]);
  
  init_grid (32);

  scalar s[];
  face vector u[];
  reset ({s,u}, undefined);
  foreach()
    s[] = 1;
  foreach_face()
    u.x[] = 1;
  
  boundary ({s,u});
  scalar * list = {s,u};

  /** the loop below is the prototype for grid/tree-common.h:refine()
      Fixes made here must also be applied to refine() */
  
  int refined, n = 0;
  do {
    refined = 0;
    tree->refined.n = 0;
    foreach_leaf()
      if (level <= 9 && sq(x - 0.5) + sq(y - 0.5) < sq(0.05)) {
	refine_cell (point, list, 0, &tree->refined);
	refined++;
	continue;
      }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (list);
      mpi_boundary_update_buffers();
      grid->tn = 0; // so that tree is not "full"
      boundary (list);
      image();
      while (balance()) {
	foreach()
	  foreach_neighbor()
            assert (s[] == 1);
	foreach_face()
	  for (int i = -2; i <= 2; i++)
	    assert ((u.x[0,i] == 1));
	image();
      }
    }
    n++;
  } while (refined);

  set_unused_pids();
  debug_mpi (qerr);

  check_restriction (s);
}
