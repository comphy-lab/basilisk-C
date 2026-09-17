/**
# Cross-spectrum transform and shell binning

`cross_shell_average()` generalizes `shell_average()` to a pair of fields.
`same` checks it reduces to the auto-spectrum when the two fields coincide.
`scaled` checks linearity in the second field. `quadrature` (same $|k|$,
$\cos$ vs $\sin$) and `offaxis` (same $|k|$, different $(k_x,k_y)$) check that
two fields on disjoint Fourier modes cross to zero *everywhere*, not just
off the naive peak bin -- the total, not one bin, is what orthogonality
constrains. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5

FILE * fp_out = NULL;

static void check (const char * tag, const char * kind, scalar a, scalar b,
                   int m, int bexp, double Eexp)
{
  int nk = nshells (m, m);
  double * pa = malloc ((size_t) m*m*sizeof(double));
  double * pb = malloc ((size_t) m*m*sizeof(double));
  int holes  = sample_scalar_plane ({a}, pa, snap_to_cell (0., 1 << ML),
                                    X0, X0 + L0, Y0, Y0 + L0, m, m);
  holes     += sample_scalar_plane ({b}, pb, snap_to_cell (0., 1 << ML),
                                    X0, X0 + L0, Y0, Y0 + L0, m, m);

  double * da = malloc (2*(size_t) m*m*sizeof(double));
  double * db = malloc (2*(size_t) m*m*sizeof(double));
  for (int i = 0; i < m*m; i++) {
    da[2*i] = pa[i], da[2*i + 1] = 0.;
    db[2*i] = pb[i], db[2*i + 1] = 0.;
  }
  fft2D_forward (da, m, m);
  fft2D_forward (db, m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  cross_shell_average (da, db, m, m, E, nk);

  double sum = 0.;
  for (int i = 0; i < nk; i++)
    sum += E[i];

  if (pid() == 0)
    fprintf (fp_out, "%s %s %d %d %d %.17g %.17g %.17g %.17g\n",
             tag, kind, m, holes, bexp, E[bexp], Eexp,
             fabs (sum - E[bexp]), sum);
  free (pa), free (pb), free (da), free (db), free (E);
}

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  int m = 1 << ML;
  init_grid (m);

  scalar a[], b[], c[], d[], e[];
  foreach() {
    a[] = 3.;                       // constant, bin 0
    b[] = cos (5.*x);               // kx = 5, bin 5
    c[] = 2.*cos (5.*x);            // same mode, scaled
    d[] = sin (5.*x);               // same |k|, quadrature phase
    e[] = cos (3.*x + 4.*y);        // |k| = 5, different (kx,ky)
  }

  if (pid() == 0) {
    fp_out = fopen ("spectra_cross.asc", "w");
    fprintf (fp_out, "# tag kind m holes bexp E Eexp leak sum\n");
  }

  check ("same",       "exact", b, b, m, 5, 0.5);   // == auto-spectrum
  check ("scaled",     "exact", b, c, m, 5, 1.0);
  check ("quadrature", "zero",  b, d, m, 5, 0.);
  check ("offaxis",    "zero",  b, e, m, 5, 0.);
  check ("constant",   "exact", a, a, m, 0, 9.0);

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py cross spectra_cross.asc 1>&2");
  }
}
