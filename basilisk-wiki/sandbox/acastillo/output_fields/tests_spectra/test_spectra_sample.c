/**
# Sampling a plane onto a lattice

`sample_scalar_plane()` on fields known everywhere, so the error *is* the
displacement between a query point and the cell centre answering it: `err_c`
from `cos(x)`, `dz` from `z`. The `face` case pins the trap `snap_to_cell()`
exists for: `z = 0` is a cell face whenever the grid size is even. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra_sample.h"

#define ML 5

FILE * fp_out = NULL;

static void report (const char * tag, scalar * list, double h, int m)
{
  int len = list_len (list);
  double * plane = malloc ((size_t) m*m*len*sizeof(double));
  int holes = sample_scalar_plane (list, plane, h,
                                   X0, X0 + L0, Y0, Y0 + L0, m, m);
  double err_c = 0., dz = 0.;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      double xi = X0 + L0*(i + 0.5)/m;
      err_c = max (err_c, fabs (plane[(i*m + j)*len]     - cos (xi)));
      dz    = max (dz,    fabs (plane[(i*m + j)*len + 2] - h));
    }
  if (pid() == 0)
    fprintf (fp_out, "%s %d %d %.17g %.17g\n", tag, m, holes, err_c, dz);
  free (plane);
}

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  init_grid (1 << ML);

  scalar c[], d[], e[];
  scalar * list = {c, d, e};

  if (pid() == 0) {
    fp_out = fopen ("spectra_sample.asc", "w");
    fprintf (fp_out, "# tag m holes err_c dz\n");
    fprintf (fp_out, "# L0 %.17g  base %d  npe %d\n", L0, 1 << ML, npe());
  }

  foreach()
    c[] = cos (x), d[] = sin (2.*y), e[] = z;

  // lattice on the grid: the query points are the cell centres
  report ("uniform", list, snap_to_cell (0., 1 << ML), 1 << ML);

  // z = 0 is a face, so the plane returned is half a cell off
  report ("face", list, 0., 1 << ML);

  // lattice coarser than the grid: displaced by half a *grid* cell
  report ("coarse", list, snap_to_cell (0., 1 << ML), 1 << (ML - 2));

#if TREE

  // refining shrinks that displacement to the finer cell size; a lattice
  // matched to the finer level removes it
  refine (sq(z) < sq(L0/8.) && level < ML + 2);
  foreach()
    c[] = cos (x), d[] = sin (2.*y), e[] = z;

  double h = snap_to_cell (0., 1 << (ML + 2));
  report ("refined", list, h, 1 << ML);
  report ("matched", list, h, 1 << (ML + 2));
#endif

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py sample spectra_sample.asc 1>&2");
  }
}
