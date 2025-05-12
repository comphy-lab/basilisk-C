#include "utils.h"
#include "output.h"
#include "check_restriction.h"

#define BGHOSTS 2

int main (int argc, char * argv[])
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : 9;

  init_grid (16);
  
  origin (-0.1, -0.5, -0.5);

  scalar s[];
  face vector u[];
  reset ({s,u}, undefined);
  foreach()
    s[] = 1;
  foreach_face()
    u.x[] = 1;
  
  boundary ({s,u});
  scalar * list = {s,u};
  
  timer t = timer_start();

#if 0  
  refine (level <= maxlevel && sq(x) + sq(y) + sq(z) < sq(0.05));
#else
  int refined;
  do {
    refined = 0;
    tree->refined.n = 0;
    foreach (serial, nowarning)
      if (is_leaf(cell) && level <= maxlevel &&
	  sq(x) + sq(y) + sq(z) < sq(0.05)) {
	refine_cell (point, list, 0, &tree->refined);
	refined++;
      }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (list);
      mpi_boundary_update_buffers();
      grid->tn = 0; // so that tree is not "full"
      boundary (list);
      foreach()
	foreach_neighbor()
	assert (s[] == 1);
      foreach_face()
	assert (u.x[] == 1);
      while (balance(0)) {
	foreach()
	  foreach_neighbor()
          assert (s[] == 1);
	foreach_face()
	  assert (u.x[] == 1);
      }
    }
  } while (refined);
#endif
  
  long nl = 0;
  foreach (serial)
    nl++;
  long nt = nl;
  mpi_all_reduce (nt, MPI_LONG, MPI_SUM);
  double elapsed = timer_elapsed(t);
  fprintf (stderr, "nl: %ld nt: %ld t: %g s: %g\n", nl, nt, elapsed, nl/elapsed);
  
  foreach()
    foreach_neighbor()
      assert (s[] == 1);
  foreach_face()
    for (int i = -2; i <= 2; i++)
      assert ((u.x[0,i] == 1));

  check_restriction (s);
  
  scalar pid[];
  foreach()
    pid[] = pid();
  output_gfs (file = "out.gfs");
}
