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
Spectra on `nz` planes evenly spaced across `[hmin, hmax]`, each snapped to a
cell centre. `format` picks the writer; `mode` applies to ASCII only, since the
HDF5 file always appends along its time axis.

Keep `[hmin, hmax]` inside the refined region and `m` at the finest grid size,
or the lattice stops landing on cell centres. Use a fixed `nz` rather than one
plane per cell, so the block shape does not change as that region grows; the
heights are in the output, so rescaling them is a post-processing choice.

`mode` defaults to append, so repeated runs in one directory accumulate blocks
-- as the profile writers do. */

void spectrum_scalar_stack (scalar * list,
                            const char * filename = "spectra.asc",
                            double hmin = 0., double hmax = 0.,
                            int nz = 1,
                            int m = N,
                            double xmin = X0, double xmax = X0 + L0,
                            double ymin = Y0, double ymax = Y0 + L0,
                            const char * mode = "a",
                            int format = SPECTRA_ASCII)
{
  int len = list_len (list), nk = nshells (m, m), holes = 0;
  double * E = malloc ((size_t) nz*len*nk*sizeof(double));
  double * z = malloc ((size_t) nz*sizeof(double));

  for (int iz = 0; iz < nz; iz++) {
    z[iz] = snap_to_cell (nz > 1 ? hmin + (hmax - hmin)*(iz + 0.5)/nz : hmin, m);
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
                            int m = N,
                            double xmin = X0, double xmax = X0 + L0,
                            double ymin = Y0, double ymax = Y0 + L0,
                            const char * mode = "a",
                            int format = SPECTRA_ASCII)
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
                            int m = N,
                            double xmin = X0, double xmax = X0 + L0,
                            double ymin = Y0, double ymax = Y0 + L0,
                            const char * mode = "a",
                            int format = SPECTRA_ASCII)
{
  spectrum_scalar_stack ({u.x, u.y, u.z}, filename, hmin, hmax, nz, m,
                         xmin, xmax, ymin, ymax, mode, format);
}

void spectrum_vector_plane (vector u,
                            const char * filename = "spectra_u.asc",
                            double h = 0.,
                            int m = N,
                            double xmin = X0, double xmax = X0 + L0,
                            double ymin = Y0, double ymax = Y0 + L0,
                            const char * mode = "a",
                            int format = SPECTRA_ASCII)
{
  spectrum_vector_stack (u, filename, h, h, 1, m,
                         xmin, xmax, ymin, ymax, mode, format);
}

