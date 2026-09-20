/**
# Sampling a plane at a coarser level than the local grid

`sample_scalar_plane_restrict()` against today's `sample_scalar_plane()` on
a grid where a band around $z = 0$ carries a refined cylinder: `region`
(`foreach_region`/`locate()`) always returns the finest cell under a query
point, so on the diagnostic lattice (one level coarser than the band) a
point inside the refined region reads an unrestricted child value, not the
coarse-cell average. `restrict` (`restriction()` + `foreach_level`) is
uniform by construction. Five fields probe different failure shapes: a
step in $x$, a linear field in $z$, a smooth periodic field with no
$z$-dependence, a $\tanh$ ramp in $z$, and an off-axis radial profile with
no symmetry to cancel the bias -- the last two checked pointwise, since
neither reduces to a single scalar keyed by $z$ alone. */

#include "grid/octree.h"
#include "utils.h"
#include "acastillo/output_fields/spectra/spectra_sample.h"
#include "acastillo/output_fields/spectra/spectra_sample_restricted.h"

#define LBASE 8
#define LBAND (LBASE + 1)
#define LCYL  (LBASE + 2)
#define NC    (1 << LBASE)

#define BANDN 24
#define BAND  (BANDN*L0/NC)

#define RCYL  (0.22*L0)
#define WRAMP (3.*L0/NC)
#define DTANH (3.*L0/NC)

// Diagnostic lattice: one level coarser than the base grid.
#define LDIAG (LBASE - 1)
#define ND    (1 << LDIAG)

#define NF 5
scalar hvar[], zlin[], smooth[], curved[], radial[];
scalar * flds;

FILE * fp_out = NULL;

double radial_profile (double x, double y)
{
  double r = sqrt (sq(x) + sq(y));
  return 0.5*(1. - tanh ((r - RCYL)/WRAMP));
}

// Analytic values at z = 0, the plane sampled below. smooth and radial
// depend on (x,y), not z, so both are checked pointwise.
double exact_hvar   (double z) { return 0.5; }
double exact_zlin   (double z) { return z;   }
double exact_smooth (double x, double y)
{
  return cos (2.*pi*x/L0)*cos (2.*pi*y/L0);
}
double exact_curved (double z) { return tanh (z/DTANH); }

/** `region`: `sample_scalar_plane()` at the diagnostic resolution --
    `locate()` returns leaves finer than the lattice inside the band. */
void region_lattice (double h, double * err, double * dmin, double * dmax,
                     int * filled, int * holes)
{
  int m1 = ND, m2 = ND, n = m1*m2;
  double plane[NF*n], dcell[n];
  for (int i = 0; i < NF*n; i++) plane[i] = nodata;
  for (int i = 0; i < n; i++) dcell[i] = nodata;

  coord box[2] = {{X0, Y0, h}, {X0 + L0, Y0 + L0, h}};
  coord ns = {m1, m2, 1};
  coord p;
  foreach_region (p, box, ns, reduction(min:plane[:NF*n])
                  reduction(min:dcell[:n])) {
    double * alias = plane;
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*ns.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*ns.y;
    alias[(i*m2 + j)*NF + 0] = hvar[];
    alias[(i*m2 + j)*NF + 1] = zlin[];
    alias[(i*m2 + j)*NF + 2] = smooth[];
    alias[(i*m2 + j)*NF + 3] = curved[];
    alias[(i*m2 + j)*NF + 4] = radial[];
    dcell[i*m2 + j] = Delta;
  }

  for (int k = 0; k < NF; k++) err[k] = 0.;
  *dmin = 1e30, *dmax = -1e30, *filled = 0, *holes = 0;
  double del = L0/m1;
  for (int i = 0; i < n; i++) {
    if (plane[i*NF] == nodata) { (*holes)++; continue; }
    double x = X0 + (i/m2 + 0.5)*del, y = Y0 + (i%m2 + 0.5)*del;
    double e[NF] = { fabs (plane[i*NF+0] - exact_hvar (h)),
                     fabs (plane[i*NF+1] - exact_zlin (h)),
                     fabs (plane[i*NF+2] - exact_smooth (x, y)),
                     fabs (plane[i*NF+3] - exact_curved (h)),
                     fabs (plane[i*NF+4] - radial_profile (x, y)) };
    for (int k = 0; k < NF; k++) if (e[k] > err[k]) err[k] = e[k];
    if (dcell[i] < *dmin) *dmin = dcell[i];
    if (dcell[i] > *dmax) *dmax = dcell[i];
    (*filled)++;
  }
}

/** `restrict`: `sample_scalar_plane_restrict()` at `LDIAG`. Uniform `Delta`
    by construction -- the fix, exercised as shipped. */
void restricted_lattice (double h, double * err, double * dmin, double * dmax,
                         int * filled, int * holes)
{
  int m1 = ND, m2 = ND, n = m1*m2;
  double plane[NF*n];

  *holes = sample_scalar_plane_restrict (flds, plane, h, LDIAG);
  double dcell = L0/ND;

  for (int k = 0; k < NF; k++) err[k] = 0.;
  *dmin = 1e30, *dmax = -1e30, *filled = 0;
  double del = L0/m1;
  for (int idx = 0; idx < n; idx++) {
    if (plane[idx*NF] == nodata) continue;
    double x = X0 + (idx/m2 + 0.5)*del, y = Y0 + (idx%m2 + 0.5)*del;
    double e[NF] = { fabs (plane[idx*NF+0] - exact_hvar (h)),
                     fabs (plane[idx*NF+1] - exact_zlin (h)),
                     fabs (plane[idx*NF+2] - exact_smooth (x, y)),
                     fabs (plane[idx*NF+3] - exact_curved (h)),
                     fabs (plane[idx*NF+4] - radial_profile (x, y)) };
    for (int k = 0; k < NF; k++) if (e[k] > err[k]) err[k] = e[k];
    if (dcell < *dmin) *dmin = dcell;
    if (dcell > *dmax) *dmax = dcell;
    (*filled)++;
  }
}

static void report (const char * tag, double * err, double dmin, double dmax,
                    int filled, int holes)
{
  if (pid() == 0)
    fprintf (fp_out, "%s %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %d\n",
             tag, holes, err[0], err[1], err[2], err[3], err[4],
             dmin, dmax, filled);
}

int main()
{
  L0 = 1.; X0 = Y0 = Z0 = -L0/2.;
  init_grid (NC);

  refine (fabs(z) < BAND && level < LBAND);
  refine (fabs(z) < BAND && sq(x) + sq(y) < sq(RCYL) && level < LCYL);

  foreach() {
    hvar[]   = (x < 0.) ? 1. : 0.;
    zlin[]   = z;
    smooth[] = cos (2.*pi*x/L0)*cos (2.*pi*y/L0);
    curved[] = tanh (z/DTANH);
    radial[] = radial_profile (x, y);
  }
  flds = {hvar, zlin, smooth, curved, radial};

  if (pid() == 0) {
    fp_out = fopen ("spectra_restrict_plane.asc", "w");
    fprintf (fp_out, "# tag holes err_hvar err_zlin err_smooth err_curved "
                    "err_radial dmin dmax filled\n");
  }

  // z = 0 is a face at even NC; snap onto the diagnostic lattice's own cell.
  double h = snap_to_cell (0., ND);

  double err1[NF], dmin1, dmax1; int filled1, holes1;
  double err2[NF], dmin2, dmax2; int filled2, holes2;
  region_lattice     (h, err1, &dmin1, &dmax1, &filled1, &holes1);
  restricted_lattice (h, err2, &dmin2, &dmax2, &filled2, &holes2);
  report ("region",   err1, dmin1, dmax1, filled1, holes1);
  report ("restrict",  err2, dmin2, dmax2, filled2, holes2);

  if (pid() == 0) {
    fclose (fp_out);
    system ("python3 ../test_spectra.py restrict_plane "
           "spectra_restrict_plane.asc 1>&2");
  }
}
