/**
# Horizontal spectra on planes

Shell-averaged spectra of a list of fields, on one or more planes of constant
$z$. For a bi-periodic domain the planes cover one period exactly, so the
transform needs no window and no detrending -- unlike the inhomogeneous
direction, where a Fourier spectrum would not be meaningful.

Sample at the finest grid size, `m = 1 << maxlevel`, on planes inside the
finest region. That is the only lattice whose points are cell centres: with a
ratio $N/m$ the query point is a cell centre only when the ratio is odd, and on
a power-of-two grid the sole odd ratio is 1. Any other choice returns one
arbitrarily chosen child cell, displaced by up to half a cell, and on an
adaptive grid that displacement varies across the plane -- position-dependent
jitter rather than a rigid shift, so it does not cancel in $|\hat c|^2$.

### Usage

```.c
  vector u[];
  spectrum_vector_stack (u, "spectra_u.h5", -0.5, 0.5, 32, 1 << MAXLEVEL,
                         X0, X0 + L0, Y0, Y0 + L0, "a", SPECTRA_HDF5);
```

see, also the tests in [tests_spectra/](../tests_spectra/): the transform and
binning in [test_spectra_modes.c](../tests_spectra/test_spectra_modes.c), the
sampling rule stated above in
[test_spectra_amr.c](../tests_spectra/test_spectra_amr.c), and the two writers
in [test_spectra_ascii.c](../tests_spectra/test_spectra_ascii.c) and
[test_spectra_hdf5.c](../tests_spectra/test_spectra_hdf5.c).
*/

#include "spectra_shell.h"
#include "spectra_output.h"

/**
Trailing parameters every orchestrator below shares: the lattice size and
domain, and the writer's mode/format. Only what precedes it -- `hmin`/`hmax`/
`nz` for a stack, `h` for a single plane -- differs between them. */

#define SPECTRA_LATTICE_ARGS \
  int m = N, \
  double xmin = X0, double xmax = X0 + L0, \
  double ymin = Y0, double ymax = Y0 + L0, \
  const char * mode = "a", \
  int format = SPECTRA_ASCII

/**
Spectra on `nz` planes evenly spaced across `[hmin, hmax]`, each snapped to a
cell centre. `format` picks the writer; `mode` applies to ASCII only, since the
HDF5 file always appends along its time axis.

Keep `[hmin, hmax]` inside the refined region and `m` at the finest grid size,
or the lattice stops landing on cell centres. `nz` may change from block to
block as that region grows: each output is a self-contained group in the
HDF5 file ([spectra_output.h](spectra_output.h)), and an ASCII block carries
its own `nz` in its header. The heights are in the output, so rescaling them
is a post-processing choice.

`mode` defaults to append, so repeated runs in one directory accumulate blocks
-- as the profile writers do. */

void spectrum_scalar_stack (scalar * list,
                            const char * filename = "spectra.asc",
                            double hmin = 0., double hmax = 0.,
                            int nz = 1,
                            SPECTRA_LATTICE_ARGS)
{
  int len = list_len (list), nk = nshells (m, m), holes = 0;
  double * E = malloc ((size_t) nz*len*nk*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));

  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (hmin + (hmax - hmin)*(iz + 0.5)/nz, m);
    holes += spectrum_plane (list, E + (size_t) iz*len*nk, z[iz],
                             xmin, xmax, ymin, ymax, m);
  }
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
  spectrum_scalar_stack (list, filename, h, h, 1, m,
                         xmin, xmax, ymin, ymax, mode, format);
}

/** Per-component spectra of a vector. The total and any horizontal/vertical
split follow by summing in post-processing. */

void spectrum_vector_stack (vector u,
                            const char * filename = "spectra_u.asc",
                            double hmin = 0., double hmax = 0.,
                            int nz = 1,
                            SPECTRA_LATTICE_ARGS)
{
  spectrum_scalar_stack ({u.x, u.y, u.z}, filename, hmin, hmax, nz, m,
                         xmin, xmax, ymin, ymax, mode, format);
}

void spectrum_vector_plane (vector u,
                            const char * filename = "spectra_u.asc",
                            double h = 0.,
                            SPECTRA_LATTICE_ARGS)
{
  spectrum_vector_stack (u, filename, h, h, 1, m,
                         xmin, xmax, ymin, ymax, mode, format);
}

/**
## Cross-spectra on planes

Transverse twin of `spectrum_scalar_stack()` for `cross_spectrum_plane()`:
`nz` planes of a single label, since the two operands need not be from a
field list the writer can name -- the same reason `cross_spectrum_foliated()`
takes a `label` below. */

void cross_spectrum_stack (scalar a, scalar b, const char * label,
                           const char * filename = "spectra_cross.asc",
                           double hmin = 0., double hmax = 0.,
                           int nz = 1,
                           SPECTRA_LATTICE_ARGS)
{
  int nk = nshells (m, m), holes = 0;
  double * E = malloc ((size_t) nz*nk*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));

  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (hmin + (hmax - hmin)*(iz + 0.5)/nz, m);
    holes += cross_spectrum_plane (a, b, E + (size_t) iz*nk, z[iz],
                                   xmin, xmax, ymin, ymax, m);
  }
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
slabs -- `sample_scalar_stack_sum()` in
[spectra_sample.h](spectra_sample.h) -- before the same FFT and shell-average
backend used above. One block, written through the same writers as
`spectrum_scalar_stack()`; `format` is unchanged, so a foliated run is told
apart from a transverse one only by its filename, which the caller chooses.

`means` follows `sample_scalar_stack_sum()`'s layout, `means[iz*len + k]` for
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
  int len = list_len (list), nk = nshells (m, m);
  size_t npt = (size_t) m*m;
  double * plane = malloc (npt*len*sizeof(double));
  double * zs = malloc ((size_t) nz*sizeof(double));
  int holes = sample_scalar_stack_sum (list, plane, means, zs,
                                       hmin, hmax, nz,
                                       xmin, xmax, ymin, ymax, m, m);
  if (holes && pid() == 0)
    fprintf (stderr,
            "spectrum_scalar_foliated: %d points unfilled over %d planes\n",
            holes, nz);
  free (zs);

  double * E = malloc ((size_t) len*nk*sizeof(double));
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
                            m, xmin, xmax, ymin, ymax, mode, format);
}

/**
## Foliated cross-spectra

`cross_spectrum_scalar_foliated()` (in
[spectra_shell.h](spectra_shell.h)) only computes; this writes the result
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
  int nk = nshells (m, m);
  double * E = malloc ((size_t) nk*sizeof(double));
  int holes = cross_spectrum_scalar_foliated (a, b, means_a, means_b, E,
                                              hmin, hmax, nz,
                                              xmin, xmax, ymin, ymax, m);
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

