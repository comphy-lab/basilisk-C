/**
# Transform and shell binning on single modes

Modes the lattice represents exactly, so the tolerances are machine
tolerances. `oblique` and `diagonal` pin the binning: |k| = 5 off-axis, and
|k| = sqrt(2), which reaches bin 1 only if the rounding in `shell_average()`
is right.

`subsampled` breaks the rule in [spectra.h](../spectra/spectra.h): halving the
lattice puts every query point on a cell face, where `locate()` picks a side
inconsistently. That jitter is not a rigid shift, so energy leaves the bin. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5

FILE * fp_out = NULL;

static void check (const char * tag, const char * kind, scalar s, int m,
                   int bexp, double Eexp, double leakexp)
{
  int nk = nshells (m, m);
  double * plane = malloc ((size_t) m*m*sizeof(double));
  int holes = sample_scalar_plane ({s}, plane, snap_to_cell (0., 1 << ML),
                                   X0, X0 + L0, Y0, Y0 + L0, m, m);
  double ms = 0.;
  for (int i = 0; i < m*m; i++)
    ms += sq (plane[i]);
  ms /= (double) m*m;

  double * data = malloc (2*(size_t) m*m*sizeof(double));
  for (int i = 0; i < m*m; i++)
    data[2*i] = plane[i], data[2*i + 1] = 0.;
  fft2D_forward (data, m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  shell_average (data, m, m, E, nk);

  double sum = 0., peak = -1.;
  int bpeak = -1;
  for (int b = 0; b < nk; b++) {
    sum += E[b];
    if (E[b] > peak)
      peak = E[b], bpeak = b;
  }
  if (pid() == 0)
    fprintf (fp_out, "%s %s %d %d %d %d %.17g %.17g %.17g %.17g %.17g\n",
             tag, kind, m, holes, bpeak, bexp, E[bexp], Eexp,
             fabs (sum - E[bexp] - leakexp), fabs (sum - ms), sum);
  free (plane), free (data), free (E);
}

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  int m = 1 << ML;
  init_grid (m);

  scalar a[], b[], d[], e[], g[];
  foreach() {
    a[] = 3.;                       // bin 0, E = 9
    b[] = cos (5.*x);               // kx = 5,        bin 5, E = 1/2
    d[] = cos (3.*x + 4.*y);        // |k| = 5,       bin 5, E = 1/2
    e[] = cos (5.*x) + 2.*sin (3.*y);  // bin 3 at 2, bin 5 at 1/2
    g[] = cos (x + y);              // |k| = sqrt(2), bin 1, E = 1/2
  }

  if (pid() == 0) {
    fp_out = fopen ("spectra_modes.asc", "w");
    fprintf (fp_out, "# tag kind m holes bpeak bexp E Eexp leak parseval sum\n");
  }

  check ("constant", "exact", a, m, 0, 9.0, 0.);
  check ("axial",    "exact", b, m, 5, 0.5, 0.);
  check ("oblique",  "exact", d, m, 5, 0.5, 0.);
  check ("twomode",  "exact", e, m, 3, 2.0, 0.5);
  check ("diagonal", "exact", g, m, 1, 0.5, 0.);

  // half the lattice: every query point on a face, so the bin bleeds
  check ("subsampled", "leaky", g, m/2, 1, 0.5, 0.);

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py modes spectra_modes.asc 1>&2");
  }
}
