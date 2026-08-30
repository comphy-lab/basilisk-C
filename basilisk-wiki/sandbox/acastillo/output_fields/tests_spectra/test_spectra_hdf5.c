/**
# The HDF5 spectrum file

Three blocks appended along the time axis, then read back in C -- which keeps
the checker on the standard library and gives it the stored values rather than
a transcription. As in [test_spectra_ascii.c](test_spectra_ascii.c) the fields
are scaled by $z$ and the mode moves between blocks, so `(nt, nz, nk)` is
checked by its contents. Only the vector file's shapes are read. */

#include "utils.h"
#include "acastillo/output_fields/spectra/spectra.h"

#define ML 5
#define NZ 4
#define NT 3

#ifdef HAVE_HDF5

// Peak bin, its energy, and the total, for one (block, plane) spectrum.
static void summarise (FILE * fp, const char * name, const double * E,
                       const double * z, int nt, int nz, int nk)
{
  for (int it = 0; it < nt; it++)
    for (int iz = 0; iz < nz; iz++) {
      const double * s = E + ((size_t) it*nz + iz)*nk;
      double sum = 0., peak = -1.;
      int bpeak = -1;
      for (int b = 0; b < nk; b++) {
        sum += s[b];
        if (s[b] > peak)
          peak = s[b], bpeak = b;
      }
      fprintf (fp, "peak %s %d %d %d %.17g %.17g %.17g\n",
               name, it, iz, bpeak, peak, sum, z[(size_t) it*nz + iz]);
    }
}

static void readback (const char * filename, const char * vecfile)
{
  FILE * fp = fopen ("spectra_hdf5.asc", "w");
  hid_t file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  double L0f;
  int mf, nzf, nkf;
  H5LTget_attribute_double (file, ".", "L0", &L0f);
  H5LTget_attribute_int (file, ".", "m",  &mf);
  H5LTget_attribute_int (file, ".", "nz", &nzf);
  H5LTget_attribute_int (file, ".", "nk", &nkf);
  fprintf (fp, "attr %.17g %d %d %d\n", L0f, mf, nzf, nkf);

  hsize_t dt[3], dz[3], de[3], dk[3];
  H5LTget_dataset_info (file, "/t",     dt, NULL, NULL);
  H5LTget_dataset_info (file, "/z",     dz, NULL, NULL);
  H5LTget_dataset_info (file, "/kphys", dk, NULL, NULL);
  H5LTget_dataset_info (file, "/E/a",   de, NULL, NULL);
  fprintf (fp, "shape t %d\n", (int) dt[0]);
  fprintf (fp, "shape z %d %d\n", (int) dz[0], (int) dz[1]);
  fprintf (fp, "shape kphys %d\n", (int) dk[0]);
  fprintf (fp, "shape E %d %d %d\n", (int) de[0], (int) de[1], (int) de[2]);

  int nt = dt[0], nz = dz[1], nk = de[2];
  double * tv = malloc (nt*sizeof(double));
  double * zv = malloc ((size_t) nt*nz*sizeof(double));
  double * kv = malloc (nk*sizeof(double));
  double * Ev = malloc ((size_t) nt*nz*nk*sizeof(double));

  H5LTread_dataset_double (file, "/t", tv);
  H5LTread_dataset_double (file, "/z", zv);
  H5LTread_dataset_double (file, "/kphys", kv);
  for (int it = 0; it < nt; it++)
    fprintf (fp, "t %d %.17g\n", it, tv[it]);

  // the bin axis is fixed, so one residual against 2 pi k / L0 says enough
  double dkmax = 0.;
  for (int b = 0; b < nk; b++)
    dkmax = max (dkmax, fabs (kv[b] - 2.*pi*b/L0));
  fprintf (fp, "kphys %.17g\n", dkmax);

  H5LTread_dataset_double (file, "/E/a", Ev);
  summarise (fp, "a", Ev, zv, nt, nz, nk);
  H5LTread_dataset_double (file, "/E/b", Ev);
  summarise (fp, "b", Ev, zv, nt, nz, nk);

  free (tv), free (zv), free (kv), free (Ev);
  H5Fclose (file);

  // the vector file goes through the same writer; check that it is there
  hid_t vfile = H5Fopen (vecfile, H5F_ACC_RDONLY, H5P_DEFAULT);
  const char * comp[3] = {"/E/v.x", "/E/v.y", "/E/v.z"};
  for (int k = 0; k < 3; k++) {
    hsize_t d[3] = {0, 0, 0};
    H5LTget_dataset_info (vfile, comp[k], d, NULL, NULL);
    fprintf (fp, "shape %s %d %d %d\n", comp[k] + 3,
             (int) d[0], (int) d[1], (int) d[2]);
  }
  H5Fclose (vfile);
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

  // three times, with a widening zone so the heights move between blocks
  for (int step = 0; step < NT; step++) {
    t = 0.1*step;
    double Lz = 0.5 + 0.25*step;
    foreach() {
      a[] = (1. + z)*cos ((5. + step)*x);   // bin 5 + step
      b[] = 3. + step + z;                  // bin 0
      v.x[] = a[], v.y[] = b[], v.z[] = cos (3.*y);
    }
    spectrum_scalar_stack ({a, b}, "spec.h5", -Lz, Lz, NZ, m,
                           X0, X0 + L0, Y0, Y0 + L0, "a", SPECTRA_HDF5);
    spectrum_vector_stack (v, "spec_u.h5", -Lz, Lz, NZ, m,
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
