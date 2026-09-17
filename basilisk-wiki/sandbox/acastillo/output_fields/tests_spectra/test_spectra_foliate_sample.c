/**
# Foliating a stack of planes

`sample_scalar_stack_sum()` on fields known everywhere, so the expected sum is
analytic. `zero_mean` supplies the exact per-slab mean of a $z$-only field, so
every slab's anomaly is 0 and the foliated sum is 0. `offset_mean` supplies a
mean of 0 for an $x$-only field, so each slab contributes the field itself and
the foliated sum is `nz` times it -- the coefficient a foliated spectrum's
normalisation has to account for. `wrong_mean` checks the anomaly is linear in
the supplied mean, not clipped or renormalised. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra_sample.h"

#define ML 5

FILE * fp_out = NULL;

static void report (const char * tag, scalar s, double * means,
                    double hmin, double hmax, int nz, int m,
                    double (* expected)(double x, double y))
{
  scalar * list = {s};
  double * plane = malloc ((size_t) m*m*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));
  int holes = sample_scalar_stack_sum (list, plane, means, z,
                                       hmin, hmax, nz,
                                       X0, X0 + L0, Y0, Y0 + L0, m, m);
  double err = 0.;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      double xi = X0 + L0*(i + 0.5)/m, yj = Y0 + L0*(j + 0.5)/m;
      err = max (err, fabs (plane[i*m + j] - expected (xi, yj)));
    }
  if (pid() == 0)
    fprintf (fp_out, "%s %d %d %d %.17g\n", tag, m, nz, holes, err);
  free (plane);
  free (z);
}

static double exp_zero (double x, double y) { return 0.; }
static double exp_cosx (double x, double y) { return cos (x); }
static double exp_cosx_m1_times4 (double x, double y) { return 4.*(cos (x) - 1.); }

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  int m = 1 << ML;
  init_grid (m);

  scalar c[], d[];
  foreach()
    c[] = z, d[] = cos (x);   // c: z-only, d: x-only (independent of z)

  int nz = 4;

  if (pid() == 0) {
    fp_out = fopen ("spectra_foliate_sample.asc", "w");
    fprintf (fp_out, "# tag m nz holes err\n");
  }

  // c's slab means are the slab heights themselves: exact per-slab mean of a
  // z-only field, so every anomaly is 0.
  {
    double means[4];
    for (int iz = 0; iz < nz; iz++)
      means[iz] = -L0/2. + L0*(iz + 0.5)/nz;
    report ("zero_mean", c, means, -L0/2., L0/2., nz, m, exp_zero);
  }

  // d does not depend on z, so a mean of 0 leaves each slab's anomaly equal
  // to d itself; the foliated sum is nz*cos(x).
  {
    double means[4] = {0., 0., 0., 0.};
    report ("offset_mean", d, means, -L0/2., L0/2., nz, m, exp_cosx);
  }

  // same field, but a deliberately wrong (nonzero) supplied mean: checks the
  // anomaly is linear in `means`, not silently corrected.
  {
    double means[4] = {1., 1., 1., 1.};
    report ("wrong_mean", d, means, -L0/2., L0/2., nz, m, exp_cosx_m1_times4);
  }

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py foliate_sample "
           "spectra_foliate_sample.asc 1>&2");
  }
}
