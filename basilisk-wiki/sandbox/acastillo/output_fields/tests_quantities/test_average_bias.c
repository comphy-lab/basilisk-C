/**
# Test for bias in average_scalar_plane
*/

#include "acastillo/output_fields/profiles/profiles.h"

int main() {
  L0 = 1.0;
  X0 = Y0 = Z0 = -0.5;
  init_grid(1 << 5); // Base level 5

  // Refine left half to level 6
  refine (x < 0 && level < 6);

  scalar s[];
  foreach()
    s[] = (x < 0) ? 1.0 : 2.0;
  boundary({s});

  // Sample over whole domain
#if dimension == 3
  coord nsamples = {100, 100, 100}; // Uniform sampling 3D
  coord region_box[2] = {{X0, Y0, Z0}, {X0 + L0, Y0 + L0, Z0 + L0}};
#else
  coord nsamples = {100, 100}; // Uniform sampling
  coord region_box[2] = {{X0, Y0}, {X0 + L0, Y0 + L0}};
#endif

  double v1 = 0;
  double v2 = 0;
  scalar * list = {s};

  average_scalar_plane(list, unity, &v1, &v2, 0., region_box, nsamples);

  fprintf(stderr, "Computed Average: %g, M2: %g\n", v1, v2);

  // Expected:
  // Volume weighted true average: (1.0 * 0.5 + 2.0 * 0.5) / 1.0 = 1.5
  // If biased by cell size (Delta):
  // Coarse (x>0): Delta_c = 2 * Delta_f. Weight 2x.
  // Avg = (1.0 * 1 + 2.0 * 2) / (1 + 2) = 5/3 = 1.666...

  if (fabs(v1 - 1.5) < 0.05 && fabs(v2 - 2.5) < 0.05) {
    fprintf(stderr, "PASS: Correct volume weighting (1.5, 2.5)\n");
    return 0; // Success
  } else if (fabs(v1 - 1.5) >= 0.05) {
    fprintf(stderr, "FAIL: Bias detected in mean (%g vs 1.5)\n", v1);
    return 1;
  } else {
    fprintf(stderr, "FAIL: Bias detected in m2 (%g vs 2.5)\n", v2);
    return 1;
  }
}
