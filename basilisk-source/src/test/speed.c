/* speed of definition of halos */

#include "utils.h"

scalar h[];

int main (int argc, char ** argv)
{
  origin (-0.5, -0.5);
  int n = 1024;
  init_grid (n);

  double R0 = 0.1;
  foreach() {
    x -= 0.5;
    h[] = exp(-(x*x + y*y)/(R0*R0));
  }
  boundary ({h});
  
  clock_t start, end0, end;
  start = end0 = clock ();
  int i;
  for (i = 0; i < 61; i++) {
    /* coarsening */
    adapt_wavelet ({h}, (double []){1e-5}, 10);
    if (i == 0)
      end0 = clock();
  }
  end = clock ();
  double cpu0 = ((double) (end0 - start))/CLOCKS_PER_SEC;
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  printf ("---- restriction + wavelet + coarsen_wavelet "
	  "+ flag_halo_cells ----\n");
  int leaves = 0, maxlevel = 0;
  foreach() { leaves++; if (level > maxlevel) maxlevel = level; }
  printf ("after coarsening: %d leaves, maximum level %d\n", 
	  leaves, maxlevel);
  printf ("initial coarsening:  %6g CPU, %.3g points.steps/s\n",
	  cpu0, n*n/cpu0);
  printf ("%4d iterations:     %6g CPU, %.3g leaves.steps/s\n",
	  i - 1, cpu - cpu0, leaves*(i - 1)/(cpu - cpu0));

  int nhalos = 0;
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l) 
      foreach_child() {
      //        printf ("%g %g %d %d\n", x, y, level, cell.neighbors);
        nhalos++;
      }

  start = clock ();
  for (i = 0; i < 10000; i++) {
    scalar * list = {h};
    int l = depth();
    boundary_iterate (level, list, 0);
    for (int i = 0; i < l; i++) {
      foreach_halo (prolongation, i)
	for (scalar s in list)
	  s.prolongation (point, s);
      boundary_iterate (level, list, i + 1);
    }
  }
  end = clock ();
  cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  printf ("---- update_halos ----\n");
  printf ("%d halo points\n", nhalos);
  printf ("%4d iterations:     %6g CPU, %.3g halos.steps/s\n",
	  i, cpu, nhalos*i/cpu);
}
