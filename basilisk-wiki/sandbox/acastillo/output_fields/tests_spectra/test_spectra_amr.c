/**
# Sampling across a refinement boundary

Same field, same lattice; only the plane moves. A level-7 slab in $z$ carries
a level-8 patch in its $x < 0$ half, so `matched` avoids the patch and every
lattice point is a cell centre, while `straddling` runs through it and half
the plane falls on faces. That jitter varies across the plane, so it does not
cancel in $|\hat c|^2$. `lmin`/`lmax` record which case is which. Tree grids
only. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5
#define KMODE 5

FILE * fp_out = NULL;

static void check (const char * tag, const char * kind, scalar s,
                   double h, int m)
{
  int lmin = 99, lmax = 0;
  foreach (reduction(min:lmin) reduction(max:lmax))
    if (fabs (z - h) < Delta/2.) {
      if (level < lmin) lmin = level;
      if (level > lmax) lmax = level;
    }

  int nk = nshells (m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  int holes = spectrum_plane ({s}, E, h, X0, X0 + L0, Y0, Y0 + L0, m);

  double sum = 0.;
  for (int b = 0; b < nk; b++)
    sum += E[b];

  // from the same sampling, so Parseval stays the transform's own identity
  double * plane = malloc ((size_t) m*m*sizeof(double));
  sample_scalar_plane ({s}, plane, h, X0, X0 + L0, Y0, Y0 + L0, m, m);
  double ms = 0.;
  for (int i = 0; i < m*m; i++)
    ms += sq (plane[i]);
  ms /= (double) m*m;

  if (pid() == 0)
    fprintf (fp_out, "%s %s %d %d %d %d %.17g %.17g %.17g\n",
             tag, kind, m, holes, lmin, lmax,
             E[KMODE], fabs (sum - E[KMODE]), fabs (sum - ms));
  free (E), free (plane);
}

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  init_grid (1 << ML);

  refine (fabs(z) < L0/6.  && level < ML + 2);              // level 7 slab
  refine (fabs(z - 0.5) < L0/32. && x < 0. && level < ML + 3); // level 8 patch

  scalar a[];
  foreach()
    a[] = cos (KMODE*x);

  if (pid() == 0) {
    fp_out = fopen ("spectra_amr.asc", "w");
    fprintf (fp_out, "# tag kind m holes lmin lmax E leak parseval\n");
  }

  check ("matched",    "exact", a, snap_to_cell (-0.5, 1 << (ML + 2)), 1 << (ML + 2));
  check ("straddling", "leaky", a, snap_to_cell ( 0.5, 1 << (ML + 3)), 1 << (ML + 2));

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py amr spectra_amr.asc 1>&2");
  }
}
