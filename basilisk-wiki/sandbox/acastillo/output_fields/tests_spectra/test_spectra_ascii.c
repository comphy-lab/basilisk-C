/**
# The ASCII spectrum block

Checks the file rather than the transform: block structure, header, and that
every value sits where it belongs. Every field is scaled by $z$, so each plane
carries a different spectrum -- transposed strides would hold the right
numbers in the wrong places. The second block differs in time, mode and
amplitude, so the append path shows too. */

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

  double hmin = -L0/8., hmax = L0/8.;
  scalar a[], b[];
  vector v[];

  // first block, written fresh: a in bin 0, b in bin 5
  t = 0.;
  foreach()
    a[] = 3. + z, b[] = (1. + z)*cos (5.*x);
  spectrum_scalar_stack ({a, b}, "spec.asc", hmin, hmax, NZ, m,
                         X0, X0 + L0, Y0, Y0 + L0, "w");

  // second block, appended: different time, amplitudes and mode
  t = 0.25;
  foreach()
    a[] = 5. + z, b[] = (2. + z)*cos (7.*x);
  spectrum_scalar_stack ({a, b}, "spec.asc", hmin, hmax, NZ, m,
                         X0, X0 + L0, Y0, Y0 + L0, "a");

  // one plane, three components, to check the column layout of a vector
  t = 0.5;
  foreach()
    v.x[] = cos (5.*x), v.y[] = 2., v.z[] = cos (3.*y);
  spectrum_vector_plane (v, "spec_u.asc", snap_to_cell (0., m), m,
                         X0, X0 + L0, Y0, Y0 + L0, "w");

  if (pid() == 0)
    system ("python3 ../test_spectra.py ascii spec.asc spec_u.asc 1>&2");
}
