// Omitting slablevel must be refused, not guessed.
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "acastillo/output_fields/profiles/profiles_slab_restrict.h"

scalar q[];

int main() {
  L0 = 1.; X0 = Y0 = Z0 = -L0/2.;
  init_grid (16);
  foreach() q[] = z;
  remove ("omitted.asc");
  profile_scalar_slab ({q}, filename = "omitted.asc", n = 16, mode = "w");
  FILE * fp = fopen ("omitted.asc", "r");
  fprintf (stderr, "file written: %s\n", fp ? "YES (BAD)" : "no (correct)");
  if (fp) fclose (fp);
}
