/**
# Shell-averaged spectra of a plane

The transform and binning behind [spectra.h](spectra.h). A plane sampled by
[spectra_sample.h](spectra_sample.h) is transformed with GSL and reduced to
`E(k)` by averaging over shells of constant $|k|$.

Bin 0 holds the squared mean, so $\sum_k E(k)$ is the mean square of the plane
and $\sum_{k>0} E(k)$ its variance -- an exact Parseval check. Bins beyond
$m/2$ are only partly populated, since they exist only in the corners of the
$k$ plane; they are kept so that the sum stays exact.
*/

#include <gsl/gsl_fft_complex.h>
#pragma autolink -lgsl -lgslcblas

#ifndef REAL // also defined by the Dimonte initial-condition headers
  #define REAL(z,i) ((z)[2*(i)])
  #define IMAG(z,i) ((z)[2*(i)+1])
#endif

#include "spectra_sample.h"

/**
Forward twin of `fft2D()` in
[initial_conditions_dimonte_fft2.h](/sandbox/acastillo/input_fields/initial_conditions_dimonte_fft2.h).
Mixed radix, so any transform length. */

static void fft2D_forward (double * data, int n0, int n1)
{
  gsl_fft_complex_wavetable * wt = gsl_fft_complex_wavetable_alloc (n1);
  gsl_fft_complex_workspace * ws = gsl_fft_complex_workspace_alloc (n1);
  for (int i = 0; i < n0; i++)
    gsl_fft_complex_forward (&REAL(data, i*n1), 1, n1, wt, ws);
  gsl_fft_complex_wavetable_free (wt);
  gsl_fft_complex_workspace_free (ws);

  wt = gsl_fft_complex_wavetable_alloc (n0);
  ws = gsl_fft_complex_workspace_alloc (n0);
  double * col = malloc (2*n0*sizeof(double));
  for (int j = 0; j < n1; j++) {
    for (int i = 0; i < n0; i++) {
      REAL(col,i) = REAL(data, i*n1 + j);
      IMAG(col,i) = IMAG(data, i*n1 + j);
    }
    gsl_fft_complex_forward (col, 1, n0, wt, ws);
    for (int i = 0; i < n0; i++) {
      REAL(data, i*n1 + j) = REAL(col,i);
      IMAG(data, i*n1 + j) = IMAG(col,i);
    }
  }
  free (col);
  gsl_fft_complex_wavetable_free (wt);
  gsl_fft_complex_workspace_free (ws);
}

// Number of shells for an m1 x m2 lattice: the corner of the k plane.
static int nshells (int m1, int m2)
{
  return (int) (sqrt (sq (m1/2.) + sq (m2/2.)) + 0.5) + 1;
}

/**
Shell-average $|F|^2/(m_1m_2)^2$ into unit-width bins in $|k|$, with $k$ in
units of $2\pi/L_0$. GSL's transform is unnormalised, hence the $1/(m_1m_2)^2$.
*/

static void shell_average (double * data, int m1, int m2, double * E, int nk)
{
  for (int b = 0; b < nk; b++)
    E[b] = 0.;
  double norm = 1./sq ((double) m1*m2);
  for (int i = 0; i < m1; i++) {
    double kx = i <= m1/2 ? i : i - m1;
    for (int j = 0; j < m2; j++) {
      double ky = j <= m2/2 ? j : j - m2;
      int b = (int) (sqrt (sq (kx) + sq (ky)) + 0.5);
      if (b < nk)
        E[b] += (sq (REAL(data, i*m2 + j)) + sq (IMAG(data, i*m2 + j)))*norm;
    }
  }
}

/**
Cross-shell-average $\mathrm{Re}[F_1 F_2^*]/(m_1m_2)^2$ into the same
unit-width $|k|$ bins as `shell_average()`. Reduces to `shell_average(data,
data, ...)` when the two transforms coincide, and to 0 on every bin where the
two fields carry power on disjoint Fourier modes -- the cross terms outside a
shared $(k_x,k_y)$ never appear in the sum.
*/

static void cross_shell_average (double * data1, double * data2,
                                 int m1, int m2, double * E, int nk)
{
  for (int b = 0; b < nk; b++)
    E[b] = 0.;
  double norm = 1./sq ((double) m1*m2);
  for (int i = 0; i < m1; i++) {
    double kx = i <= m1/2 ? i : i - m1;
    for (int j = 0; j < m2; j++) {
      double ky = j <= m2/2 ? j : j - m2;
      int b = (int) (sqrt (sq (kx) + sq (ky)) + 0.5);
      if (b < nk)
        E[b] += (REAL(data1, i*m2 + j)*REAL(data2, i*m2 + j) +
                 IMAG(data1, i*m2 + j)*IMAG(data2, i*m2 + j))*norm;
    }
  }
}

/**
Sample the plane at $z = h$ and fill `E` with one shell-averaged spectrum per
field, laid out as `E[is*nk + b]` for field `is` in bin `b`. The caller
allocates `len*nk` doubles, `nk` from `nshells()`. Returns the number of
lattice points no rank owned, 0 for a complete plane.

`m` must be the finest grid size, or the plane is not sampled on cell centres.
Only rank 0 transforms, since the sampled plane is replicated everywhere.
*/

int spectrum_plane (scalar * list, double * E, double h,
                    double xmin, double xmax, double ymin, double ymax, int m)
{
  int len = list_len (list), nk = nshells (m, m);
  size_t npt = (size_t) m*m;
  double * plane = malloc (npt*len*sizeof(double));
  int holes = sample_scalar_plane (list, plane, h,
                                   xmin, xmax, ymin, ymax, m, m);
  if (pid() == 0) {
    double * data = malloc (2*npt*sizeof(double));
    for (int is = 0; is < len; is++) {
      for (size_t i = 0; i < npt; i++) {
        REAL(data,i) = plane[i*len + is];
        IMAG(data,i) = 0.;
      }
      fft2D_forward (data, m, m);
      shell_average (data, m, m, E + (size_t) is*nk, nk);
    }
    free (data);
  }
  free (plane);
  return holes/len;   // per lattice point, not per stored value
}

/**
Cross-spectrum twin of `spectrum_plane()`: sample two fields at $z = h$ and
fill `E` with their shell-averaged $\mathrm{Re}[\hat a \hat b^*]$. Reduces to
an auto-spectrum when `b` and `a` are the same field. `holes` is reported from
`a`'s sample only -- it depends on domain coverage, not field values, so both
samples return the same count (see `sample_scalar_plane()`).
*/

int cross_spectrum_plane (scalar a, scalar b, double * E, double h,
                          double xmin, double xmax, double ymin, double ymax,
                          int m)
{
  int nk = nshells (m, m);
  size_t npt = (size_t) m*m;
  double * pa = malloc (npt*sizeof(double));
  double * pb = malloc (npt*sizeof(double));
  int holes = sample_scalar_plane ({a}, pa, h, xmin, xmax, ymin, ymax, m, m);
  sample_scalar_plane ({b}, pb, h, xmin, xmax, ymin, ymax, m, m);
  if (pid() == 0) {
    double * da = malloc (2*npt*sizeof(double));
    double * db = malloc (2*npt*sizeof(double));
    for (size_t i = 0; i < npt; i++) {
      REAL(da,i) = pa[i], IMAG(da,i) = 0.;
      REAL(db,i) = pb[i], IMAG(db,i) = 0.;
    }
    fft2D_forward (da, m, m);
    fft2D_forward (db, m, m);
    cross_shell_average (da, db, m, m, E, nk);
    free (da);
    free (db);
  }
  free (pa);
  free (pb);
  return holes;
}

/**
Foliated twin of `cross_spectrum_plane()`: sum each field's anomaly over `nz`
slabs with `sample_scalar_stack_sum()` (see
[spectra_sample.h](spectra_sample.h)), then cross-shell-average the two
resulting planes. `means_a`/`means_b` follow that function's `means`
layout, one per field. `holes` is the sum over both fields and all slabs.
*/

int cross_spectrum_scalar_foliated (scalar a, scalar b,
                                    const double * means_a,
                                    const double * means_b,
                                    double * E,
                                    double hmin, double hmax, int nz,
                                    double xmin, double xmax,
                                    double ymin, double ymax, int m)
{
  int nk = nshells (m, m);
  size_t npt = (size_t) m*m;
  double * pa = malloc (npt*sizeof(double));
  double * pb = malloc (npt*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));
  int holes = sample_scalar_stack_sum ({a}, pa, means_a, z, hmin, hmax, nz,
                                       xmin, xmax, ymin, ymax, m, m);
  holes += sample_scalar_stack_sum ({b}, pb, means_b, z, hmin, hmax, nz,
                                    xmin, xmax, ymin, ymax, m, m);
  if (pid() == 0) {
    double * da = malloc (2*npt*sizeof(double));
    double * db = malloc (2*npt*sizeof(double));
    for (size_t i = 0; i < npt; i++) {
      REAL(da,i) = pa[i], IMAG(da,i) = 0.;
      REAL(db,i) = pb[i], IMAG(db,i) = 0.;
    }
    fft2D_forward (da, m, m);
    fft2D_forward (db, m, m);
    cross_shell_average (da, db, m, m, E, nk);
    free (da);
    free (db);
  }
  free (pa);
  free (pb);
  free (z);
  return holes;
}
