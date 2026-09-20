/**
# Round-robin FFT jobs against the serial reference

`spectrum_scalar_stack()` and `cross_spectrum_stack()`
([spectra_restricted.h](../spectra/spectra_restricted.h)) sample every plane
collectively, then hand the independent FFT + shell-average jobs that follow
out round-robin across ranks (`job % npe() == pid()`) instead of running them
all on `pid() == 0`. Since every rank already holds the same sampled plane
(via the sampling reduction) and only the small `E` array is summed back
together at the end, the result must be exactly independent of *which* rank
does each job -- so `E`, read back from what the real functions write, must
be bit-identical to a serial reference computed with `pid() == 0` doing every
job, run under the same `npe()`. Disagreeing here means a job was
double-counted or dropped, not a numerical difference. */

#include "grid/octree.h"
#include "utils.h"
#include "acastillo/output_fields/spectra/spectra_restricted.h"

#define ML 6
#define NZ 8

FILE * fp_out = NULL;

/** Serial reference for `spectrum_scalar_stack()`: every job on
    `pid() == 0`, inlined so no file I/O is needed for this half. */
static void serial_stack (scalar * list, double * E, double * z,
                          double hmin, double hmax, int nz, int slablevel)
{
  int m = 1 << slablevel, len = list_len (list), nk = nshells (m, m);
  double * plane = malloc ((size_t) m*m*len*sizeof(double));
  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (nz > 1 ? hmin + (hmax - hmin)*(iz + 0.5)/nz : hmin, m);
    sample_scalar_plane_restrict (list, plane, z[iz], slablevel);
    if (pid() == 0) {
      fft2D_plan plan = fft2D_plan_new (m);
      double * data = malloc (2*(size_t) m*m*sizeof(double));
      for (int is = 0; is < len; is++) {
        for (size_t i = 0; i < (size_t) m*m; i++) {
          REAL(data,i) = plane[i*len + is];
          IMAG(data,i) = 0.;
        }
        fft2D_forward (data, plan);
        shell_average (data, m, m, E + ((size_t) iz*len + is)*nk, nk);
      }
      free (data);
      fft2D_plan_free (plan);
    }
  }
  free (plane);
}

/** Serial reference for `cross_spectrum_stack()`: one `cross_spectrum_plane()`
    call per slab, `pid() == 0` only inside it. */
static void cross_stack_serial (scalar a, scalar b, double * E, double * z,
                                double hmin, double hmax, int nz,
                                int slablevel)
{
  int m = 1 << slablevel, nk = nshells (m, m);
  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (nz > 1 ? hmin + (hmax - hmin)*(iz + 0.5)/nz : hmin, m);
    cross_spectrum_plane (a, b, E + (size_t) iz*nk, z[iz], slablevel);
  }
}

/** Read back `spectrum_scalar_stack()`'s own ASCII block: `nz*nk` rows,
    `len` energy columns starting at column 4. Skips the two header lines. */
static void read_stack_block (const char * filename, double * E, int nz,
                              int nk, int len)
{
  FILE * fp = fopen (filename, "r");
  if (fp == NULL) { perror (filename); exit (1); }
  char line[4096];
  fgets (line, sizeof line, fp);   // "# Spectrum: ..."
  fgets (line, sizeof line, fp);   // "# [0]iz ..."
  for (int iz = 0; iz < nz; iz++)
    for (int b = 0; b < nk; b++) {
      int iz_r, b_r; double z_r, k_r;
      if (fscanf (fp, "%d %lf %d %lf", &iz_r, &z_r, &b_r, &k_r) != 4) {
        fprintf (stderr, "%s: short read at iz=%d b=%d\n", filename, iz, b);
        exit (1);
      }
      for (int is = 0; is < len; is++)
        if (fscanf (fp, "%lf", &E[((size_t) iz*len + is)*nk + b]) != 1) {
          fprintf (stderr, "%s: short read at iz=%d b=%d is=%d\n",
                  filename, iz, b, is);
          exit (1);
        }
    }
  fclose (fp);
}

/** Read back `cross_spectrum_stack()`'s own ASCII block: `nz*nk` rows, one
    energy column at column 4. */
static void read_cross_stack_block (const char * filename, double * E,
                                    int nz, int nk)
{
  FILE * fp = fopen (filename, "r");
  if (fp == NULL) { perror (filename); exit (1); }
  char line[4096];
  fgets (line, sizeof line, fp);   // "# Cross-spectrum(...): ..."
  fgets (line, sizeof line, fp);   // "# [0]iz ..."
  for (int iz = 0; iz < nz; iz++)
    for (int b = 0; b < nk; b++) {
      int iz_r, b_r; double z_r, k_r;
      if (fscanf (fp, "%d %lf %d %lf %lf", &iz_r, &z_r, &b_r, &k_r,
                 &E[(size_t) iz*nk + b]) != 5) {
        fprintf (stderr, "%s: short read at iz=%d b=%d\n", filename, iz, b);
        exit (1);
      }
    }
  fclose (fp);
}

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  init_grid (1 << ML);

  scalar a[], b[], d[];
  foreach() {
    a[] = 3.;                          // bin 0
    b[] = cos (5.*x);                  // kx = 5, bin 5
    d[] = cos (3.*x + 4.*y);           // |k| = 5, bin 5
  }
  scalar * list = {a, b, d};
  int len = list_len (list);

  int slablevel = ML;
  int m = 1 << slablevel, nk = nshells (m, m);
  double hmin = -0.3*L0, hmax = 0.3*L0;
  int nz = NZ;

  double * Eref = malloc ((size_t) nz*len*nk*sizeof(double));
  double * zref = malloc ((size_t) nz*sizeof(double));
  double * Err  = malloc ((size_t) nz*len*nk*sizeof(double));

  serial_stack (list, Eref, zref, hmin, hmax, nz, slablevel);
  spectrum_scalar_stack (list, "spectra_restrict_rr_stack.asc", hmin, hmax,
                        nz, slablevel);
  if (pid() == 0)
    read_stack_block ("spectra_restrict_rr_stack.asc", Err, nz, nk, len);

  double maxdiff = 0.;
  if (pid() == 0)
    for (size_t i = 0; i < (size_t) nz*len*nk; i++) {
      double diff = fabs (Err[i] - Eref[i]);
      if (diff > maxdiff) maxdiff = diff;
    }

  double * Ecref = malloc ((size_t) nz*nk*sizeof(double));
  double * zcref = malloc ((size_t) nz*sizeof(double));
  double * Ecrr  = malloc ((size_t) nz*nk*sizeof(double));

  cross_stack_serial (b, d, Ecref, zcref, hmin, hmax, nz, slablevel);
  cross_spectrum_stack (b, d, "bd", "spectra_restrict_rr_cross.asc",
                       hmin, hmax, nz, slablevel);
  if (pid() == 0)
    read_cross_stack_block ("spectra_restrict_rr_cross.asc", Ecrr, nz, nk);

  double maxcdiff = 0.;
  if (pid() == 0)
    for (size_t i = 0; i < (size_t) nz*nk; i++) {
      double diff = fabs (Ecrr[i] - Ecref[i]);
      if (diff > maxcdiff) maxcdiff = diff;
    }

  if (pid() == 0) {
    fp_out = fopen ("spectra_restrict_rr.asc", "w");
    fprintf (fp_out, "# tag npe nz len m maxdiff\n");
    fprintf (fp_out, "scalar_stack %d %d %d %d %.17g\n",
             npe(), nz, len, m, maxdiff);
    fprintf (fp_out, "cross_stack %d %d %d %d %.17g\n",
             npe(), nz, 2, m, maxcdiff);
    fclose (fp_out);
    system ("python3 ../test_spectra.py restrict_rr spectra_restrict_rr.asc 1>&2");
  }

  free (Eref), free (zref), free (Err);
  free (Ecref), free (zcref), free (Ecrr);
}
