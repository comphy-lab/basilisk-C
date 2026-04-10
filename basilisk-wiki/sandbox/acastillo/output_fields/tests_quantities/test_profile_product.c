/**
# Test profile_product_foreach_region

We test the profile product function by computing the product of two linear fields.
In 2D: s1 = y, s2 = y. Product should be y^2.
In 3D: s1 = z, s2 = z. Product should be z^2.
*/

#include "grid/octree.h"
#include "acastillo/output_fields/profiles/profiles.h"

int main() {
  L0 = 1.0;
  X0 = Y0 = Z0 = -0.5;
  init_grid(16);

  scalar s1[], s2[];

  foreach() {
#if dimension == 2
    s1[] = y;
    s2[] = y;
#else
    s1[] = z;
    s2[] = z;
#endif
  }

  scalar * list1 = {s1};
  scalar * list2 = {s2};

  // We expect the profile to be exactly z^2 (or y^2).
  // The profile is computed at specific heights.
  // The averaging is done over x (and y in 3D).
  // Since fields only depend on the profile coordinate (y or z),
  // the average should be EXACTLY the value at that height.

  char * filename = "profiles_prod.asc";

  // Compute profile using default arguments for geometry
  double del = L0 / N;
  double hmin = -L0/2. + del/2.;
  double hmax =  L0/2. - del/2.;
  profile_foreach_region(list1, unity, "profiles.asc", hmin=hmin, hmax=hmax);
  profile_product_foreach_region(list1, list2, unity, filename, hmin=hmin, hmax=hmax);

  // Verification
  if (pid() == 0) {
    FILE * fp = fopen(filename, "r");
    if (fp == NULL) {
      perror(filename);
      return 1;
    }

    char line[256];
    // Skip header lines (start with '#')
    while (fgets(line, sizeof(line), fp) && line[0] == '#');

    double max_err = 0;
    int count = 0;
    int iprof;
    double h, delta, val;

    // Format: iprof, h, delta, product
    do {
      if (sscanf(line, "%d %lf %lf %lf", &iprof, &h, &delta, &val) == 4) {
        double expected = h * h;
        double err = fabs(val - expected);
        if (err > max_err) max_err = err;
        count++;
      }
    } while (fgets(line, sizeof(line), fp));
    fclose(fp);

    fprintf(stderr, "Processed %d points. Max Error: %g\n", count, max_err);

    if (count == 0 || max_err > 1e-9) {
      fprintf(stderr, "FAIL: Profile product mismatch.\n");
      return 1;
    } else {
      fprintf(stderr, "PASS: Profile product correct.\n");
      return 0;
    }
  }

  return 0;
}
