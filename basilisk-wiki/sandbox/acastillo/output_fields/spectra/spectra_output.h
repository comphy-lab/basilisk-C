/**
# Writing a spectrum block

Two writers over the same array: `E[(iz*len + is)*nk + b]` for plane `iz`,
field `is`, bin `b`, with `z` holding the `nz` plane heights. Both run on rank
0 only, which is where [spectra_shell.h](spectra_shell.h) leaves the data.
*/

#define SPECTRA_ASCII 0
#define SPECTRA_HDF5  1

// HDF5 is optional; the ASCII writer is always available.
#if __has_include(<hdf5.h>)
  #ifndef HAVE_HDF5
    #define HAVE_HDF5 1
  #endif
  #pragma autolink -lhdf5 -lhdf5_hl
  #include <hdf5.h>
  #include <hdf5_hl.h>
#endif

/**
Write one block: two header lines, then `nz*nk` rows ordered by plane and then
by bin, then a blank line. Rank 0 only. `E` holds `E[(iz*len + is)*nk + b]` for
plane `iz`, field `is`, bin `b`, and `z` the `nz` plane heights. */

void write_spectrum_block_ascii (const char * filename, const char * mode,
                           scalar * list, double * E, double * z,
                           int nz, int nk, int m, double hmin, double hmax)
{
  if (pid() != 0)
    return;

  int len = list_len (list);
  FILE * fp = fopen (filename, mode);
  if (fp == NULL) { perror (filename); exit (1); }

  fprintf (fp, "# Spectrum: t = %.10g, L0 = %g, m = %d, nz = %d, nk = %d,"
           " hmin = %g, hmax = %g\n", t, L0, m, nz, nk, hmin, hmax);
  fprintf (fp, "# [0]iz [1]z [2]k [3]kphys");
  int c = 4;
  for (scalar s in list)
    fprintf (fp, " [%d]E(%s)", c++, s.name);
  fputc ('\n', fp);

  for (int iz = 0; iz < nz; iz++)
    for (int b = 0; b < nk; b++) {
      fprintf (fp, "%-4d %15.8e %-6d %15.8e", iz, z[iz], b, 2.*pi*b/L0);
      for (int is = 0; is < len; is++)
        fprintf (fp, " %24.15e", E[((size_t) iz*len + is)*nk + b]);
      fputc ('\n', fp);
    }

  fputc ('\n', fp);
  fflush (fp);
  fclose (fp);
}

/**
## A single cross-spectrum block

`cross_spectrum_scalar_foliated()` and `cross_spectrum_plane()` return one
shell-averaged array with no field list to draw a name from -- the two
operands are arbitrary scalars, not necessarily `u.x`/`u.y`/`u.z` -- so the
caller supplies a `label` string instead. One block: `nk` rows, no `z`/`iz`
column since there is exactly one plane (or one foliated sum) per call.
*/

void write_cross_spectrum_block_ascii (const char * filename,
                                       const char * mode,
                                       const char * label,
                                       double * E, int nk,
                                       double hmin, double hmax)
{
  if (pid() != 0)
    return;

  FILE * fp = fopen (filename, mode);
  if (fp == NULL) { perror (filename); exit (1); }

  fprintf (fp, "# Cross-spectrum(%s): t = %.10g, L0 = %g, nk = %d,"
           " hmin = %g, hmax = %g\n", label, t, L0, nk, hmin, hmax);
  fprintf (fp, "# [0]k [1]kphys [2]E(%s)\n", label);
  for (int b = 0; b < nk; b++)
    fprintf (fp, "%-6d %15.8e %24.15e\n", b, 2.*pi*b/L0, E[b]);

  fputc ('\n', fp);
  fflush (fp);
  fclose (fp);
}

/**
## A cross-spectrum stack

Transverse twin of `write_spectrum_block_ascii()` for a single label spread
over `nz` planes: `E[iz*nk + b]`, one `z` per plane -- the cross-spectrum
counterpart of `spectrum_scalar_stack()`.
*/

void write_cross_spectrum_stack_ascii (const char * filename, const char * mode,
                                       const char * label,
                                       double * E, double * z,
                                       int nz, int nk, double hmin, double hmax)
{
  if (pid() != 0)
    return;

  FILE * fp = fopen (filename, mode);
  if (fp == NULL) { perror (filename); exit (1); }

  fprintf (fp, "# Cross-spectrum(%s): t = %.10g, L0 = %g, nz = %d, nk = %d,"
           " hmin = %g, hmax = %g\n", label, t, L0, nz, nk, hmin, hmax);
  fprintf (fp, "# [0]iz [1]z [2]k [3]kphys [4]E(%s)\n", label);
  for (int iz = 0; iz < nz; iz++)
    for (int b = 0; b < nk; b++)
      fprintf (fp, "%-4d %15.8e %-6d %15.8e %24.15e\n",
               iz, z[iz], b, 2.*pi*b/L0, E[(size_t) iz*nk + b]);

  fputc ('\n', fp);
  fflush (fp);
  fclose (fp);
}

/**
## HDF5

One group per output, holding everything that output produced:

    /                      attr: format = 3
    /t0000.8500/           attrs: t, hmin, hmax, L0, nz, nk (and m)
        k        (nk)      bin index
        kphys    (nk)      $2\pi k/L_0$
        z        (nz)      plane heights, snapped to cell centres
        <field>  (nz, nk)  one per field, named after it

Every dataset is created at its final shape and written once, contiguous and
unfiltered. Nothing is extendible and nothing is appended: `H5DOappend()`
reads the trailing extents from the dataset rather than from the buffer, so
a block whose plane count had changed ran off the end of the caller's array,
and it aborts inside HDF5 1.12 even where the shapes do agree.

`nz` and `nk` are therefore free to change between outputs, and a group
carries its own axes instead of sharing them -- a few kB of repetition per
output, against no invariant a later block can break. A reader opens one
group and has everything.

The group is named for its time, zero-padded so that lexicographic order is
time order; the exact value is the `t` attribute. Re-running a time replaces
its group, since a restart recomputes it.

One file per run rather than one per output, unlike the field writers: a
block is a few hundred kB, so the whole run fits in one file. Serial HDF5,
since rank 0 already holds the whole array -- `create_hdf5_file()` in
[output_hdf5_helpers.h](/sandbox/acastillo/output_fields/output_hdf5_helpers.h)
is collective and would hang if called from one rank.
*/

#ifdef HAVE_HDF5

#define SPECTRA_HDF5_FORMAT 3

// One create, one write, at the shape the caller already has. `ncols` of 0
// means rank 1.
static void write_fixed_dataset (hid_t loc, const char * name,
                                 const void * data, hid_t type,
                                 hsize_t nrows, hsize_t ncols)
{
  if (!nrows)
    return;
  hsize_t dims[2] = {nrows, ncols};
  hid_t sp = H5Screate_simple (ncols ? 2 : 1, dims, NULL);
  hid_t dset = H5Dcreate2 (loc, name, type, sp, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
  if (dset < 0)
    fprintf (stderr, "spectra: cannot create %s\n", name);
  else {
    if (H5Dwrite (dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
      fprintf (stderr, "spectra: cannot write %s\n", name);
    H5Dclose (dset);
  }
  H5Sclose (sp);
}

// The run's file, created on the first block. Probed with stdio rather than
// H5E_BEGIN_TRY, which qcc cannot parse. A file from an older layout is
// refused rather than grown a second set of groups beside its own.
static hid_t spectra_open_file (const char * filename, const char * who)
{
  bool created = true;
  FILE * probe = fopen (filename, "r");
  if (probe) {
    fclose (probe);
    created = false;
  }
  hid_t file = created ?
    H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) :
    H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    fprintf (stderr, "%s: cannot open %s\n", who, filename);
    return -1;
  }

  int format = SPECTRA_HDF5_FORMAT;
  if (created)
    H5LTset_attribute_int (file, ".", "format", &format, 1);
  else {
    int fmt = -1;
    if (H5LTget_attribute_int (file, ".", "format", &fmt) < 0 ||
        fmt != SPECTRA_HDF5_FORMAT) {
      fprintf (stderr, "%s: %s is layout %d, this writer is %d -- delete it, "
               "or keep it and write elsewhere\n", who, filename, fmt, format);
      H5Fclose (file);
      return -1;
    }
  }
  return file;
}

// This output's group, with the axes and metadata every writer stores.
static hid_t spectra_block_group (hid_t file, double hmin, double hmax,
                                  int nz, int nk)
{
  char name[64];
  snprintf (name, sizeof(name), "t%09.4f", t);
  if (H5Lexists (file, name, H5P_DEFAULT) > 0)
    H5Ldelete (file, name, H5P_DEFAULT);
  hid_t g = H5Gcreate2 (file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (g < 0) {
    fprintf (stderr, "spectra: cannot create group %s\n", name);
    return -1;
  }

  double tt = t;
  H5LTset_attribute_double (g, ".", "t", &tt, 1);
  H5LTset_attribute_double (g, ".", "hmin", &hmin, 1);
  H5LTset_attribute_double (g, ".", "hmax", &hmax, 1);
  H5LTset_attribute_double (g, ".", "L0", &L0, 1);
  H5LTset_attribute_int (g, ".", "nz", &nz, 1);
  H5LTset_attribute_int (g, ".", "nk", &nk, 1);

  int * kidx = malloc (nk*sizeof(int));
  double * kphys = malloc (nk*sizeof(double));
  for (int b = 0; b < nk; b++) {
    kidx[b] = b;
    kphys[b] = 2.*pi*b/L0;
  }
  write_fixed_dataset (g, "k", kidx, H5T_NATIVE_INT, nk, 0);
  write_fixed_dataset (g, "kphys", kphys, H5T_NATIVE_DOUBLE, nk, 0);
  free (kidx);
  free (kphys);

  return g;
}

void write_spectrum_block_hdf5 (const char * filename,
                                scalar * list, double * E, double * z,
                                int nz, int nk, int m,
                                double hmin, double hmax)
{
  if (pid() != 0)
    return;

  hid_t file = spectra_open_file (filename, "write_spectrum_block_hdf5");
  if (file < 0)
    return;
  hid_t g = spectra_block_group (file, hmin, hmax, nz, nk);
  if (g < 0) {
    H5Fclose (file);
    return;
  }
  H5LTset_attribute_int (g, ".", "m", &m, 1);
  write_fixed_dataset (g, "z", z, H5T_NATIVE_DOUBLE, nz, 0);

  // gather each field, which is strided in E, into a contiguous nz x nk block
  int len = list_len (list);
  double * slice = malloc ((size_t) nz*nk*sizeof(double));
  int is = 0;
  for (scalar sc in list) {
    for (int iz = 0; iz < nz; iz++)
      for (int b = 0; b < nk; b++)
        slice[(size_t) iz*nk + b] = E[((size_t) iz*len + is)*nk + b];
    write_fixed_dataset (g, sc.name, slice, H5T_NATIVE_DOUBLE, nz, nk);
    is++;
  }
  free (slice);

  H5Gclose (g);
  H5Fclose (file);
}

/** Transverse twin of `write_spectrum_block_hdf5()` for a single label over
`nz` planes: one dataset named after the label instead of one per field,
beside the same `z` axis. */

void write_cross_spectrum_stack_hdf5 (const char * filename,
                                      const char * label,
                                      double * E, double * z,
                                      int nz, int nk,
                                      double hmin, double hmax)
{
  if (pid() != 0)
    return;

  hid_t file = spectra_open_file (filename,
                                  "write_cross_spectrum_stack_hdf5");
  if (file < 0)
    return;
  hid_t g = spectra_block_group (file, hmin, hmax, nz, nk);
  if (g < 0) {
    H5Fclose (file);
    return;
  }

  write_fixed_dataset (g, "z", z, H5T_NATIVE_DOUBLE, nz, 0);
  write_fixed_dataset (g, label, E, H5T_NATIVE_DOUBLE, nz, nk);

  H5Gclose (g);
  H5Fclose (file);
}

/** Same layout again for a folded block, which is a single row: `z` holds
the midpoint of the folded range, the height `spectrum_scalar_foliated()`
stores for the same reason. */

void write_cross_spectrum_block_hdf5 (const char * filename,
                                      const char * label,
                                      double * E, int nk,
                                      double hmin, double hmax)
{
  if (pid() != 0)
    return;

  hid_t file = spectra_open_file (filename,
                                  "write_cross_spectrum_block_hdf5");
  if (file < 0)
    return;
  hid_t g = spectra_block_group (file, hmin, hmax, 1, nk);
  if (g < 0) {
    H5Fclose (file);
    return;
  }

  double zc = 0.5*(hmin + hmax);
  write_fixed_dataset (g, "z", &zc, H5T_NATIVE_DOUBLE, 1, 0);
  write_fixed_dataset (g, label, E, H5T_NATIVE_DOUBLE, 1, nk);

  H5Gclose (g);
  H5Fclose (file);
}

#endif // HAVE_HDF5
