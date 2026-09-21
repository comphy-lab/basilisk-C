/**
# The HDF5 spectrum file

Three blocks written to one file, then read back in C -- which keeps the
checker on the standard library and gives it the stored values rather than a
transcription. As in [test_spectra_ascii.c](test_spectra_ascii.c) the fields
are scaled by $z$ and the mode moves between blocks, so the stored array is
checked by its contents.

Each output is its own group, so the groups are enumerated rather than
assumed: the test reads whatever the writer left, in the order HDF5 lists it
by name, and checks that order is time order -- which is what the zero
padding in the group name buys.

The plane count changes between blocks (`NZS`), as it does for a caller
sampling a zone that grows. A fixed `nz` here would not exercise it, and did
not: the shape was wrong on disk for two layouts running.

Only the vector file's shapes are read. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5
#define NT 3

// planes per block: a different count each time, including one that shrinks
static const int NZS[NT] = {4, 2, 7};

#ifdef HAVE_HDF5

// Peak bin, its energy, and the total, for each plane of one spectrum.
static void summarise (FILE * fp, const char * name, const double * E,
                       const double * z, int it, int nz, int nk)
{
  for (int iz = 0; iz < nz; iz++) {
    const double * s = E + (size_t) iz*nk;
    double sum = 0., peak = -1.;
    int bpeak = -1;
    for (int b = 0; b < nk; b++) {
      sum += s[b];
      if (s[b] > peak)
        peak = s[b], bpeak = b;
    }
    fprintf (fp, "peak %s %d %d %d %.17g %.17g %.17g\n",
             name, it, iz, bpeak, peak, sum, z[iz]);
  }
}

// Every group of one file, in the order HDF5 lists them by name.
static void dump_file (FILE * fp, const char * filename,
                       const char ** fields, int nfields, bool values)
{
  hid_t file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    fprintf (fp, "missing %s\n", filename);
    return;
  }

  if (values) {
    int format = -1;
    H5LTget_attribute_int (file, ".", "format", &format);
    fprintf (fp, "format %d\n", format);
  }

  // H5Gget_num_objs() rather than H5Gget_info(), whose H5G_info_t qcc
  // cannot parse. Root holds nothing but the per-output groups.
  hsize_t ngroups = 0;
  H5Gget_num_objs (file, &ngroups);
  for (hsize_t i = 0; i < ngroups; i++) {
    char gname[64];
    H5Lget_name_by_idx (file, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                        gname, sizeof(gname), H5P_DEFAULT);
    hid_t g = H5Gopen2 (file, gname, H5P_DEFAULT);
    int it = i;

    hsize_t dz[2] = {0, 0}, de[2] = {0, 0};
    H5LTget_dataset_info (g, "z", dz, NULL, NULL);
    double * zv = malloc (dz[0]*sizeof(double));
    H5LTread_dataset_double (g, "z", zv);

    if (values) {
      double tv, hminv, hmaxv;
      int nzv, nkv, mv;
      H5LTget_attribute_double (g, ".", "t", &tv);
      H5LTget_attribute_double (g, ".", "hmin", &hminv);
      H5LTget_attribute_double (g, ".", "hmax", &hmaxv);
      H5LTget_attribute_int (g, ".", "nz", &nzv);
      H5LTget_attribute_int (g, ".", "nk", &nkv);
      H5LTget_attribute_int (g, ".", "m", &mv);
      fprintf (fp, "block %d %s %d %d %d %.17g %.17g %.17g\n",
               it, gname, nzv, nkv, mv, tv, hminv, hmaxv);
      fprintf (fp, "shape %d z %d\n", it, (int) dz[0]);

      // the bin axis is stored per group, so check each one
      hsize_t dk[2] = {0, 0};
      H5LTget_dataset_info (g, "kphys", dk, NULL, NULL);
      double * kv = malloc (dk[0]*sizeof(double));
      H5LTread_dataset_double (g, "kphys", kv);
      double dkmax = 0.;
      for (int b = 0; b < (int) dk[0]; b++)
        dkmax = max (dkmax, fabs (kv[b] - 2.*pi*b/L0));
      fprintf (fp, "kphys %d %d %.17g\n", it, (int) dk[0], dkmax);
      free (kv);
    }

    for (int k = 0; k < nfields; k++) {
      H5LTget_dataset_info (g, fields[k], de, NULL, NULL);
      fprintf (fp, "shape %d %s %d %d\n", it, fields[k],
               (int) de[0], (int) de[1]);
      if (values) {
        double * Ev = malloc ((size_t) de[0]*de[1]*sizeof(double));
        H5LTread_dataset_double (g, fields[k], Ev);
        summarise (fp, fields[k], Ev, zv, it, (int) de[0], (int) de[1]);
        free (Ev);
      }
    }

    free (zv);
    H5Gclose (g);
  }
  H5Fclose (file);
}

static void readback (const char * filename, const char * vecfile)
{
  FILE * fp = fopen ("spectra_hdf5.asc", "w");
  const char * scalars[2] = {"a", "b"};
  const char * comps[3] = {"v.x", "v.y", "v.z"};
  dump_file (fp, filename, scalars, 2, true);
  // the vector file goes through the same writer; check that it is there
  dump_file (fp, vecfile, comps, 3, false);
  fclose (fp);
}

#endif // HAVE_HDF5

int main()
{
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  int m = 1 << ML;
  init_grid (m);

  scalar a[], b[];
  vector v[];

  // three times, with a widening zone so the heights move between blocks,
  // and a different plane count each time
  for (int step = 0; step < NT; step++) {
    t = 0.1*step;
    double Lz = 0.5 + 0.25*step;
    foreach() {
      a[] = (1. + z)*cos ((5. + step)*x);   // bin 5 + step
      b[] = 3. + step + z;                  // bin 0
      v.x[] = a[], v.y[] = b[], v.z[] = cos (3.*y);
    }
    spectrum_scalar_stack ({a, b}, "spec.h5", -Lz, Lz, NZS[step], m,
                           X0, X0 + L0, Y0, Y0 + L0, "a", SPECTRA_HDF5);
    spectrum_vector_stack (v, "spec_u.h5", -Lz, Lz, NZS[step], m,
                           X0, X0 + L0, Y0, Y0 + L0, "a", SPECTRA_HDF5);
  }

#ifdef HAVE_HDF5
  if (pid() == 0) {
    readback ("spec.h5", "spec_u.h5");
    system ("python3 ../test_spectra.py hdf5 spectra_hdf5.asc 1>&2");
  }
#else
  if (pid() == 0)
    fprintf (stderr, "HDF5 not available\n");
#endif
}
