/**
# Horizontal spectra on planes, restriction-sampled

Level-based twin of [spectra.h](spectra.h). Shell-averaged spectra of a list
of fields, on one or more planes of constant $z$. For a bi-periodic domain
the planes cover one period exactly, so the transform needs no window and no
detrending -- unlike the inhomogeneous direction, where a Fourier spectrum
would not be meaningful.

Sample on the `slablevel` lattice ($m = 2^{\text{slablevel}}$) via
`restriction()` + `foreach_level(slablevel)`
([spectra_sample_restricted.h](spectra_sample_restricted.h)): every cell
contributes its own restricted value at exactly that level, so the lattice is
uniform whether or not the interface refines past `slablevel` -- unlike
[spectra.h](spectra.h)'s `locate()`-based sampling, which required the
lattice to be the finest grid level or it would land off cell centres.

### Usage

```.c
  vector u[];
  spectrum_vector_stack (u, "spectra_u.h5", -0.5, 0.5, 32, MAXLEVEL,
                         "a", SPECTRA_HDF5);
```

see, also the tests in [tests_spectra/](../tests_spectra/): the transform and
binning in [test_spectra_modes.c](../tests_spectra/test_spectra_modes.c), the
sampling rule stated above in
[test_spectra_amr.c](../tests_spectra/test_spectra_amr.c), and the two writers
in [test_spectra_ascii.c](../tests_spectra/test_spectra_ascii.c) and
[test_spectra_hdf5.c](../tests_spectra/test_spectra_hdf5.c).
*/

#include "spectra_shell_restricted.h"
#include "spectra_output.h"

/**
Trailing parameters every orchestrator below shares: the diagnostic level and
the writer's mode/format. Only what precedes it -- `hmin`/`hmax`/`nz` for a
stack, `h` for a single plane -- differs between them. */

#define SPECTRA_LATTICE_ARGS \
  int slablevel, \
  const char * mode = "a", \
  int format = SPECTRA_ASCII

/**
Spectra on `nz` planes evenly spaced across `[hmin, hmax]`, each snapped to a
cell centre. `format` picks the writer; `mode` applies to ASCII only, since the
HDF5 file always appends along its time axis.

`slablevel` sets the lattice, `m = 2^slablevel` -- `restriction()` +
`foreach_level(slablevel)` (`sample_scalar_plane_restrict()` in
[spectra_sample_restricted.h](spectra_sample_restricted.h)), correct whether
`slablevel` is the finest grid level or coarser, unlike
[spectra.h](spectra.h)'s `locate()`-based sampling which required it to be
the finest. `nz` may change from block to block as the refined region grows:
each output is a self-contained group in the HDF5 file
([spectra_output.h](spectra_output.h)), and an ASCII block carries its own
`nz` in its header. The heights are in the output, so rescaling them is a
post-processing choice.

`mode` defaults to append, so repeated runs in one directory accumulate blocks
-- as the profile writers do.

Sampling stays collective; the `nz*len` FFT jobs after it are independent,
so they run round-robin (index `iz*len + is`) instead of serially on
`pid() == 0`, and `E` is reduced across ranks at the end. */

void spectrum_scalar_stack (scalar * list,
                            const char * filename = "spectra.asc",
                            double hmin = 0., double hmax = 0.,
                            int nz = 1,
                            SPECTRA_LATTICE_ARGS)
{
  int m = 1 << slablevel, len = list_len (list), nk = nshells (m, m), holes = 0;
  size_t npt = (size_t) m*m;
  double * E = malloc ((size_t) nz*len*nk*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));
  for (size_t b = 0; b < (size_t) nz*len*nk; b++)
    E[b] = 0.;

  fft2D_plan plan = fft2D_plan_new (m);
  double * plane = malloc (npt*len*sizeof(double));
  double * data = malloc (2*npt*sizeof(double));
  int job = 0;
  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (hmin + (hmax - hmin)*(iz + 0.5)/nz, m);
    holes += sample_scalar_plane_restrict (list, plane, z[iz], slablevel);

    for (int is = 0; is < len; is++, job++) {
      if (job % npe() != pid()) continue;
      for (size_t i = 0; i < npt; i++) {
        REAL(data,i) = plane[i*len + is];
        IMAG(data,i) = 0.;
      }
      fft2D_forward (data, plan);
      shell_average (data, m, m, E + ((size_t) iz*len + is)*nk, nk);
    }
  }
  free (data);
  free (plane);
  fft2D_plan_free (plan);
  mpi_all_reduce_array (E, MPI_DOUBLE, MPI_SUM, nz*len*nk);

  if (holes && pid() == 0)
    fprintf (stderr, "spectrum_scalar_stack: %d points unfilled over %d planes\n",
             holes, nz);

#ifdef HAVE_HDF5
  if (format == SPECTRA_HDF5)
    write_spectrum_block_hdf5 (filename, list, E, z, nz, nk, m, hmin, hmax);
  else
#endif
    write_spectrum_block_ascii (filename, mode, list, E, z, nz, nk, m,
                                hmin, hmax);
  free (E);
  free (z);
}

/** A single plane, as one row block of the above. */

void spectrum_scalar_plane (scalar * list,
                            const char * filename = "spectra.asc",
                            double h = 0.,
                            SPECTRA_LATTICE_ARGS)
{
  spectrum_scalar_stack (list, filename, h, h, 1, slablevel, mode, format);
}

/** Per-component spectra of a vector. The total and any horizontal/vertical
split follow by summing in post-processing. */

void spectrum_vector_stack (vector u,
                            const char * filename = "spectra_u.asc",
                            double hmin = 0., double hmax = 0.,
                            int nz = 1,
                            SPECTRA_LATTICE_ARGS)
{
  spectrum_scalar_stack ({u.x, u.y, u.z}, filename, hmin, hmax, nz,
                         slablevel, mode, format);
}

void spectrum_vector_plane (vector u,
                            const char * filename = "spectra_u.asc",
                            double h = 0.,
                            SPECTRA_LATTICE_ARGS)
{
  spectrum_vector_stack (u, filename, h, h, 1, slablevel, mode, format);
}

/**
## Cross-spectra on planes

Transverse twin of `spectrum_scalar_stack()` for `cross_spectrum_plane()`:
`nz` planes of a single label, since the two operands need not be from a
field list the writer can name -- the same reason `cross_spectrum_foliated()`
takes a `label` below.

Only 2 fields per plane, so the `nz` planes themselves are handed out
round-robin (`iz % npe() == pid()`) instead of an `is` axis; sampling stays
collective, `E` is reduced across ranks at the end. */

void cross_spectrum_stack (scalar a, scalar b, const char * label,
                           const char * filename = "spectra_cross.asc",
                           double hmin = 0., double hmax = 0.,
                           int nz = 1,
                           SPECTRA_LATTICE_ARGS)
{
  int m = 1 << slablevel, nk = nshells (m, m), holes = 0;
  size_t npt = (size_t) m*m;
  double * E = malloc ((size_t) nz*nk*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));
  for (size_t i = 0; i < (size_t) nz*nk; i++)
    E[i] = 0.;

  fft2D_plan plan = fft2D_plan_new (m);
  double * pa = malloc (npt*sizeof(double));
  double * pb = malloc (npt*sizeof(double));
  double * da = malloc (2*npt*sizeof(double));
  double * db = malloc (2*npt*sizeof(double));
  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (hmin + (hmax - hmin)*(iz + 0.5)/nz, m);
    holes += sample_scalar_plane_restrict ({a}, pa, z[iz], slablevel);
    sample_scalar_plane_restrict ({b}, pb, z[iz], slablevel);

    if (iz % npe() != pid()) continue;
    for (size_t i = 0; i < npt; i++) {
      REAL(da,i) = pa[i], IMAG(da,i) = 0.;
      REAL(db,i) = pb[i], IMAG(db,i) = 0.;
    }
    fft2D_forward (da, plan);
    fft2D_forward (db, plan);
    cross_shell_average (da, db, m, m, E + (size_t) iz*nk, nk);
  }
  free (pa);
  free (pb);
  free (da);
  free (db);
  fft2D_plan_free (plan);
  mpi_all_reduce_array (E, MPI_DOUBLE, MPI_SUM, nz*nk);

  if (holes && pid() == 0)
    fprintf (stderr, "cross_spectrum_stack(%s): %d points unfilled over %d planes\n",
             label, holes, nz);

#ifdef HAVE_HDF5
  if (format == SPECTRA_HDF5)
    write_cross_spectrum_stack_hdf5 (filename, label, E, z, nz, nk, hmin, hmax);
  else
#endif
    write_cross_spectrum_stack_ascii (filename, mode, label, E, z, nz, nk,
                                      hmin, hmax);
  free (E);
  free (z);
}

/**
## Foliated spectra

The $z$-integrated spectrum of Poujade & Peybernes (2010) and Soulard (2024):
each field's anomaly to a caller-supplied per-slab mean is summed over `nz`
slabs -- `sample_scalar_stack_sum_restrict()` in
[spectra_sample_restricted.h](spectra_sample_restricted.h) -- before the same
FFT and shell-average backend used above. One block, written through the same
writers as `spectrum_scalar_stack()`; `format` is unchanged, so a foliated
run is told apart from a transverse one only by its filename, which the
caller chooses.

`means` follows `sample_scalar_stack_sum_restrict()`'s layout, `means[iz*len + k]` for
field `k` at slab `iz` -- the per-slab mean the caller computed, by whatever
profile reduction it uses. The stored height is the midpoint of `[hmin,
hmax]`, a representative location for the zone the spectrum was folded over;
`hmin`/`hmax` themselves are in the block header, as they already are for a
transverse stack.
*/

void spectrum_scalar_foliated (scalar * list, const double * means,
                               const char * filename = "spectra_foliated.asc",
                               double hmin = 0., double hmax = 0.,
                               int nz = 1,
                               SPECTRA_LATTICE_ARGS)
{
  int m = 1 << slablevel, len = list_len (list), nk = nshells (m, m);
  size_t npt = (size_t) m*m;
  double * plane = malloc (npt*len*sizeof(double));
  double * zs = malloc ((size_t) nz*sizeof(double));
  int holes = sample_scalar_stack_sum_restrict (list, plane, means, zs,
                                                hmin, hmax, nz, slablevel);
  if (holes && pid() == 0)
    fprintf (stderr,
            "spectrum_scalar_foliated: %d points unfilled over %d planes\n",
            holes, nz);
  free (zs);

  double * E = malloc ((size_t) len*nk*sizeof(double));
  if (pid() == 0) {
    fft2D_plan plan = fft2D_plan_new (m);
    double * data = malloc (2*npt*sizeof(double));
    for (int is = 0; is < len; is++) {
      for (size_t i = 0; i < npt; i++) {
        REAL(data,i) = plane[i*len + is];
        IMAG(data,i) = 0.;
      }
      fft2D_forward (data, plan);
      shell_average (data, m, m, E + (size_t) is*nk, nk);
    }
    free (data);
    fft2D_plan_free (plan);
  }
  free (plane);

  double zc[1] = {0.5*(hmin + hmax)};
#ifdef HAVE_HDF5
  if (format == SPECTRA_HDF5)
    write_spectrum_block_hdf5 (filename, list, E, zc, 1, nk, m, hmin, hmax);
  else
#endif
    write_spectrum_block_ascii (filename, mode, list, E, zc, 1, nk, m,
                                hmin, hmax);
  free (E);
}

/** Per-component foliated spectra of a vector, as `spectrum_vector_stack()`
is to `spectrum_scalar_stack()`. */

void spectrum_vector_foliated (vector u, const double * means,
                               const char * filename = "spectra_u_foliated.asc",
                               double hmin = 0., double hmax = 0.,
                               int nz = 1,
                               SPECTRA_LATTICE_ARGS)
{
  spectrum_scalar_foliated ({u.x, u.y, u.z}, means, filename, hmin, hmax, nz,
                            slablevel, mode, format);
}

/**
## Foliated cross-spectra

`cross_spectrum_scalar_foliated()` (in
[spectra_shell_restricted.h](spectra_shell_restricted.h)) only computes; this writes the result
through the same two writers as the single-field orchestrators above, keyed
by a caller-supplied `label` since a cross term has no field of its own to
name it after -- e.g. `label = "cz"` for $(f, u_z)$, or `"st_x"` for
$(F^\sigma_x, u_x)$, one call per component so the total follows by summing
in post-processing.
*/

void cross_spectrum_foliated (scalar a, scalar b,
                              const double * means_a, const double * means_b,
                              const char * label,
                              const char * filename = "spectra_cross_foliated.asc",
                              double hmin = 0., double hmax = 0.,
                              int nz = 1,
                              SPECTRA_LATTICE_ARGS)
{
  int m = 1 << slablevel, nk = nshells (m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  int holes = cross_spectrum_scalar_foliated (a, b, means_a, means_b, E,
                                              hmin, hmax, nz, slablevel);
  if (holes && pid() == 0)
    fprintf (stderr,
            "cross_spectrum_foliated(%s): %d points unfilled over %d planes\n",
            label, holes, nz);

#ifdef HAVE_HDF5
  if (format == SPECTRA_HDF5)
    write_cross_spectrum_block_hdf5 (filename, label, E, nk, hmin, hmax);
  else
#endif
    write_cross_spectrum_block_ascii (filename, mode, label, E, nk,
                                      hmin, hmax);
  free (E);
}

#undef SPECTRA_LATTICE_ARGS

