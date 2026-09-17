/**
# Foliated spectra

`spectrum_scalar_foliated()` feeds `sample_scalar_stack_sum()`'s output into
the same FFT/shell-average backend as `spectrum_scalar_stack()`, then writes
one block through the same writer -- so what needs checking here is the
foliation, not the transform or the file format (both already covered by
[test_spectra_modes.c](test_spectra_modes.c) and
[test_spectra_ascii.c](test_spectra_ascii.c)).

`exact_mean` supplies each slab's true mean for a $z$-only field, so the
foliated sum is exactly 0 and every bin is 0. `zero_mean` folds a $z$-independent
mode with a mean of 0, so the foliated sum is `nz` times the mode and
$E$(bin 5) $= (nz)^2/2$. `wrong_mean` reuses that mode with a deliberately
nonzero mean, picking up an extra constant that must land in bin 0 --
checking the foliation is linear in the supplied mean end to end, through the
writer. `cross_self` checks `cross_spectrum_scalar_foliated()` against the
same mode crossed with itself: it must reproduce `zero_mean`'s $E$(bin 5). */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5
#define NZ 4

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  int m = 1 << ML;
  init_grid (m);

  double hmin = -L0/2., hmax = L0/2.;
  scalar a[], b[], c[];
  foreach() {
    a[] = z;                // z-only: exact_mean's anomaly is 0 everywhere
    b[] = cos (5.*x);       // z-independent mode, kx = 5
    c[] = cos (5.*x);       // same mode, crossed against a deliberately wrong mean
  }

  double means_exact[NZ];
  for (int iz = 0; iz < NZ; iz++)
    means_exact[iz] = hmin + (hmax - hmin)*(iz + 0.5)/NZ;
  double means_zero[NZ] = {0., 0., 0., 0.};
  double means_wrong[NZ] = {1., 1., 1., 1.};

  t = 0.;
  const char * fn = "spectra_foliate.asc";
  spectrum_scalar_foliated ({a}, means_exact, fn, hmin, hmax, NZ, m,
                            X0, X0 + L0, Y0, Y0 + L0, "w");
  spectrum_scalar_foliated ({b}, means_zero, fn, hmin, hmax, NZ, m,
                            X0, X0 + L0, Y0, Y0 + L0, "a");
  spectrum_scalar_foliated ({c}, means_wrong, fn, hmin, hmax, NZ, m,
                            X0, X0 + L0, Y0, Y0 + L0, "a");

  int nk = nshells (m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  int holes = cross_spectrum_scalar_foliated (b, b, means_zero, means_zero, E,
                                              hmin, hmax, NZ,
                                              X0, X0 + L0, Y0, Y0 + L0, m);
  if (pid() == 0) {
    FILE * fp = fopen ("spectra_foliate_cross.asc", "w");
    fprintf (fp, "# holes bexp E Eexp\n");
    fprintf (fp, "%d %d %.17g %.17g\n", holes, 5, E[5], 0.5*sq ((double) NZ));
    fclose (fp);
  }
  free (E);

  if (pid() == 0)
    system ("python3 ../test_spectra.py foliate "
           "spectra_foliate.asc spectra_foliate_cross.asc 1>&2");
}
