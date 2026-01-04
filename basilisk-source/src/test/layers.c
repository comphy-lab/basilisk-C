/**
# Layers */

#define LAYERS 1
// fixme: does not work on quadtrees
#include "grid/multigrid.h"

int main()
{
  init_grid (1);
  nl = 5;
  scalar a[], h = new scalar[nl];
  reset ({a,h}, 2.);
  foreach (serial)
    foreach_layer()
      fprintf (stderr, "%g %g\n", h[], a[]);
  foreach()
    foreach_layer() {
      a[] = 5;
      h[] = point.l;
    }
  foreach (serial)
    foreach_layer()
      fprintf (stderr, "%g %g\n", h[], a[]);
  scalar hleft = new scalar[nl];
  reset ({hleft}, 33);
  foreach()
    foreach_layer()
      hleft[] = h[-1];
  foreach (serial)
    foreach_layer()
      fprintf (stderr, "%g\n", hleft[]);
  int l = 1;
  foreach_layer() {
    foreach()
      h[] = l;
    l *= 2;
  }
  foreach (serial)
    foreach_layer()
      fprintf (stderr, "%g\n", h[]);  
  delete ({h, hleft});
}
