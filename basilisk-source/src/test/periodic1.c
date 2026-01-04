#define BGHOSTS 2

int main() {
  foreach_dimension()
    periodic (right);
  init_grid (8);
  size (4);
  origin (-L0/2., -L0/2.);
#if 0
  foreach_cell_all()
    fprintf (stderr, "%g %g %d %p\n", x, y, level, cellp);
#else
  scalar a[];
  foreach()
    a[] = x + y;

  output_cells (stdout);
  foreach()
    foreach_neighbor(2)
      fprintf (stderr, "%g %g %g\n", x, y, a[]);
#endif
}
