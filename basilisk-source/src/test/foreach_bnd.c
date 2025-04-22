/* Check that foreach_boundary() works when two layers of ghost cells
   are used. */

#define BGHOSTS 2

int main() {
  init_grid (16);
  foreach_boundary (left)
    fprintf (stderr, "%g %g\n", x, y);
}
