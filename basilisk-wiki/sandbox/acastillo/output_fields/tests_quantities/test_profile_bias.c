/**
# Test for bias in profile_foreach_region
*/

#define LEVEL 5
#include "acastillo/output_fields/profiles/profiles.h"

int main() {
  L0 = 1.0;
  X0 = Y0 = Z0 = -0.5;
  init_grid(1 << LEVEL);

  // Refine left half (x < 0)
#if TREE  
  refine (x < 0 && level < 6);
#endif

  scalar s[];
  foreach()
    s[] = (x < 0) ? 1.0 : 2.0;
  boundary({s});

  scalar * list = {s};

  // Profile along Y (averaging over X at each Y)
  // Correct volume-weighted average is 1.5.
#if dimension == 3
  double hmin = Z0 + 0.1, hmax = Z0 + 0.9;
#else
  double hmin = Y0 + 0.1, hmax = Y0 + 0.9;
#endif
  profile_foreach_region(list, unity, "profiles_bias.asc", xmin = -L0/2., xmax = L0/2., hmin = hmin, hmax = hmax);

  // Check results
  if (pid() == 0) {
    FILE * fp = fopen("profiles_bias.asc", "r");
    char line[256];
    // Skip header lines (start with '#')
    while (fgets(line, sizeof(line), fp) && line[0] == '#');

    double max_err = 0;
    int iprof;
    double y_val, delta, s_mean, s_sq;
    // Format: iprof, y, delta, mean, mean_sq
    do {
      if (sscanf(line, "%d %lf %lf %lf %lf", &iprof, &y_val, &delta, &s_mean, &s_sq) == 5) {
        double err_mean = fabs(s_mean - 1.5);
        if (err_mean > max_err) max_err = err_mean;

        // Expected mean_sq = (1^2 + 2^2)/2 = 2.5
        double err_sq = fabs(s_sq - 2.5);
        if (err_sq > max_err) max_err = err_sq;
      }
    } while (fgets(line, sizeof(line), fp));
    fclose(fp);

    fprintf(stderr, "Max Error: %g\n", max_err);
    if (max_err > 0.05) {
      fprintf(stderr, "FAIL: Bias detected in profile.\n");
      return 1;
    } else {
      fprintf(stderr, "PASS: Profile unbiased.\n");
      return 0;
    }
  }
  return 0;
}
