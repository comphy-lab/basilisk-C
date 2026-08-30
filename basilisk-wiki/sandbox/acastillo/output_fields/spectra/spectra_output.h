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
## HDF5

The same block, appended along an unlimited time axis of a single file:

    /t          (nt)            simulation time
    /hmin /hmax (nt)            requested bounds, which move with the zone
    /k /kphys   (nk)            bin index and $2\pi k/L_0$
    /z          (nt, nz)        plane heights, snapped to cell centres
    /E/<field>  (nt, nz, nk)

One file per run rather than one per output, unlike the field writers: a block
is a few hundred kB, so the whole run fits in one file and the analysis can
slice it by time, height or wavenumber. Serial HDF5, since rank 0 already holds
the whole array -- `create_hdf5_file()` in
[output_hdf5_helpers.h](/sandbox/acastillo/output_fields/output_hdf5_helpers.h)
is collective and would hang if called from one rank.
*/

#ifdef HAVE_HDF5

// Create an extendible, chunked, compressed dataset of shape (0, d1, d2).
// H5LT only makes fixed-size datasets, so creation stays manual; appending
// afterwards is one H5DOappend call.
static hid_t open_or_create_extendible (hid_t file, const char * name,
                                        hid_t type, int rank,
                                        hsize_t d1, hsize_t d2,
                                        int compression_level)
{
  if (H5Lexists (file, name, H5P_DEFAULT) > 0)
    return H5Dopen2 (file, name, H5P_DEFAULT);

  hsize_t dims[3] = {0, d1, d2}, maxdims[3] = {H5S_UNLIMITED, d1, d2};
  hsize_t chunk[3] = {1, d1, d2};
  hid_t sp = H5Screate_simple (rank, dims, maxdims);
  hid_t pl = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk (pl, rank, chunk);
  H5Pset_shuffle (pl);
  H5Pset_deflate (pl, compression_level);
  hid_t dset = H5Dcreate2 (file, name, type, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
  H5Pclose (pl);
  H5Sclose (sp);
  return dset;
}

// One time slice onto the unlimited axis 0.
static void append_slice (hid_t file, const char * name, const void * data,
                          hid_t type, int rank, hsize_t d1, hsize_t d2,
                          int compression_level)
{
  hid_t dset = open_or_create_extendible (file, name, type, rank, d1, d2,
                                          compression_level);
  H5DOappend (dset, H5P_DEFAULT, 0, 1, type, data);
  H5Dclose (dset);
}

void write_spectrum_block_hdf5 (const char * filename,
                                scalar * list, double * E, double * z,
                                int nz, int nk, int m,
                                double hmin, double hmax,
                                int compression_level = 6)
{
  if (pid() != 0)
    return;

  int len = list_len (list);
  // Probe with stdio rather than H5E_BEGIN_TRY, which qcc cannot parse.
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
    fprintf (stderr, "write_spectrum_block_hdf5: cannot open %s\n", filename);
    return;
  }

  if (created) {                       // fixed axes and metadata, written once
    H5LTset_attribute_double (file, ".", "L0", &L0, 1);
    H5LTset_attribute_int (file, ".", "m",  &m,  1);
    H5LTset_attribute_int (file, ".", "nz", &nz, 1);
    H5LTset_attribute_int (file, ".", "nk", &nk, 1);

    int * kidx = malloc (nk*sizeof(int));
    double * kphys = malloc (nk*sizeof(double));
    for (int b = 0; b < nk; b++) {
      kidx[b] = b;
      kphys[b] = 2.*pi*b/L0;
    }
    hsize_t dk = nk;
    H5LTmake_dataset_int (file, "/k", 1, &dk, kidx);
    H5LTmake_dataset_double (file, "/kphys", 1, &dk, kphys);
    free (kidx);
    free (kphys);

    H5Gclose (H5Gcreate2 (file, "/E", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  }

  double tt = t;
  append_slice (file, "/t",    &tt,   H5T_NATIVE_DOUBLE, 1, 1, 1, compression_level);
  append_slice (file, "/hmin", &hmin, H5T_NATIVE_DOUBLE, 1, 1, 1, compression_level);
  append_slice (file, "/hmax", &hmax, H5T_NATIVE_DOUBLE, 1, 1, 1, compression_level);
  append_slice (file, "/z",    z,     H5T_NATIVE_DOUBLE, 2, nz, 1, compression_level);

  // gather each field, which is strided in E, into a contiguous nz x nk slice
  double * slice = malloc ((size_t) nz*nk*sizeof(double));
  int is = 0;
  for (scalar sc in list) {
    for (int iz = 0; iz < nz; iz++)
      for (int b = 0; b < nk; b++)
        slice[(size_t) iz*nk + b] = E[((size_t) iz*len + is)*nk + b];
    char name[128];
    snprintf (name, sizeof(name), "/E/%s", sc.name);
    append_slice (file, name, slice, H5T_NATIVE_DOUBLE, 3, nz, nk,
                  compression_level);
    is++;
  }
  free (slice);
  H5Fclose (file);
}

#endif // HAVE_HDF5
