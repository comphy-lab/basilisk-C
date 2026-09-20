/**
# Sampling strategies on a mixed-resolution grid

Fields with closed-form horizontal averages, sampled every way the outputs do.
Reports the error against the analytic answer and the spread of `Delta` over
the cells that answered:

1. `foreach()`, unweighted       -- `profile_*_slab` today
2. `foreach()`, `dv()`-weighted  -- the histograms
3. `foreach_region` lattice      -- the spectra
4. `restriction` + `foreach_level` -- one level, uniform `Delta`

Grid: LBASE everywhere, a band about $z=0$ at LBAND, a cylinder through the
band at LCYL. A slab in the band spans two jumps and a curved inner boundary.
The cylinder axis is $z$, so slab means there are height-independent.
*/

#include "grid/octree.h"
#include "utils.h"

#define LBASE 6
#define LBAND (LBASE + 1)
#define LCYL  (LBASE + 2)
#define NC    (1 << LBASE)

// Band half-width, in base cells.
#define BANDN 6
#define BAND  (BANDN*L0/NC)

// Cylinder radius (axis along z), and the width of the smooth radial ramp.
#define RCYL  (0.22*L0)
#define WRAMP (3.*L0/NC)

// tanh thickness, in base cells
#define DTANH (3.*L0/NC)

scalar hvar[], zlin[], smooth[], curved[], radial[];
scalar * flds;

/** `radial` varies exactly where the finest cells are. Smooth rather than a
sharp disc: a staircased boundary would swamp the sampling error measured
here. */

double radial_profile (double x, double y)
{
  double r = sqrt (sq(x) + sq(y));
  return 0.5*(1. - tanh ((r - RCYL)/WRAMP));
}

// Analytic slab means.
double exact_hvar   (double z) { return 0.5; }
double analytic_zlin(double z) { return z;   }
double exact_smooth (double z) { return 0.0; }
double exact_curved (double z) { return tanh (z/DTANH); }

/** Slab mean of `radial`, by quadrature. */
double exact_radial (double z)
{
  static double cached = -1.;
  if (cached < 0.) {
    int nq = 4096;
    double h = L0/nq, s = 0.;
    for (int i = 0; i < nq; i++)
      for (int j = 0; j < nq; j++)
        s += radial_profile (X0 + (i + 0.5)*h, Y0 + (j + 0.5)*h);
    cached = s/((double) nq*nq);
  }
  return cached;
}

#define NF 5
typedef struct { double err[NF]; double dmin, dmax; int n; } Result;

/** Cell-centre z of slab j. */
double slab_z (int j) { return Z0 + (j + 0.5)*L0/NC; }

/** 1: leaves, unweighted. 2: leaves, dv()-weighted. */
void leaves (int weighted, Result * r)
{
  int n = NC;
  double num[NF*n], den[n], dmn[n], dmx[n];
  for (int i = 0; i < NF*n; i++) num[i] = 0.;
  for (int i = 0; i < n; i++) { den[i] = 0.; dmn[i] = 1e30; dmx[i] = -1e30; }

  foreach (reduction(+:num[:NF*n]) reduction(+:den[:n])
           reduction(min:dmn[:n]) reduction(max:dmx[:n])) {
    int j = (int)((z - Z0)/(L0/NC));
    if (j < 0) j = 0; if (j >= n) j = n - 1;
    double w = weighted ? dv() : 1.;
    den[j] += w;
    num[NF*j + 0] += w*hvar[];
    num[NF*j + 1] += w*zlin[];
    num[NF*j + 2] += w*smooth[];
    num[NF*j + 3] += w*curved[];
    num[NF*j + 4] += w*radial[];
    if (Delta < dmn[j]) dmn[j] = Delta;
    if (Delta > dmx[j]) dmx[j] = Delta;
  }

  for (int k = 0; k < NF; k++) r->err[k] = 0.;
  r->dmin = 1e30; r->dmax = -1e30; r->n = 0;
  for (int j = 0; j < n; j++) {
    if (den[j] <= 0.) continue;
    double zc = slab_z (j);
    double e[NF] = { fabs (num[NF*j+0]/den[j] - exact_hvar (zc)),
                     fabs (num[NF*j+1]/den[j] - analytic_zlin (zc)),
                     fabs (num[NF*j+2]/den[j] - exact_smooth (zc)),
                     fabs (num[NF*j+3]/den[j] - exact_curved (zc)),
                     fabs (num[NF*j+4]/den[j] - exact_radial (zc)) };
    for (int k = 0; k < NF; k++) if (e[k] > r->err[k]) r->err[k] = e[k];
    if (fabs (zc) < BAND) {
      if (dmn[j] < r->dmin) r->dmin = dmn[j];
      if (dmx[j] > r->dmax) r->dmax = dmx[j];
    }
    r->n++;
  }
}

/** 3: foreach_region lattice, NC x NC per slab. */
void lattice (Result * r)
{
  for (int k = 0; k < NF; k++) r->err[k] = 0.;
  r->dmin = 1e30; r->dmax = -1e30; r->n = 0;

  for (int j = 0; j < NC; j++) {
    double zc = slab_z (j);
    coord box[2] = {{X0, Y0, zc}, {X0 + L0, Y0 + L0, zc}};
    coord ns = {NC, NC, 1};
    double s[NF] = {0.,0.,0.,0.,0.}, cnt = 0., dmn = 1e30, dmx = -1e30;
    coord p;
    foreach_region (p, box, ns, reduction(+:s[:NF]) reduction(+:cnt)
                    reduction(min:dmn) reduction(max:dmx)) {
      s[0] += hvar[]; s[1] += zlin[]; s[2] += smooth[]; s[3] += curved[];
      s[4] += radial[];
      cnt++;
      if (Delta < dmn) dmn = Delta;
      if (Delta > dmx) dmx = Delta;
    }
    if (cnt <= 0.) continue;
    double e[NF] = { fabs (s[0]/cnt - exact_hvar (zc)),
                     fabs (s[1]/cnt - analytic_zlin (zc)),
                     fabs (s[2]/cnt - exact_smooth (zc)),
                     fabs (s[3]/cnt - exact_curved (zc)),
                     fabs (s[4]/cnt - exact_radial (zc)) };
    for (int k = 0; k < NF; k++) if (e[k] > r->err[k]) r->err[k] = e[k];
    if (fabs (zc) < BAND) {
      if (dmn < r->dmin) r->dmin = dmn;
      if (dmx > r->dmax) r->dmax = dmx;
    }
    r->n++;
  }
}

/** 4: restriction, then foreach_level. Uniform Delta by design. */
void restricted (Result * r)
{
  restriction (flds);

  int n = NC;
  double num[NF*n], den[n], dmn[n], dmx[n];
  for (int i = 0; i < NF*n; i++) num[i] = 0.;
  for (int i = 0; i < n; i++) { den[i] = 0.; dmn[i] = 1e30; dmx[i] = -1e30; }

  foreach_level (LBASE, reduction(+:num[:NF*n]) reduction(+:den[:n])
                 reduction(min:dmn[:n]) reduction(max:dmx[:n])) {
    int j = (int)((z - Z0)/(L0/NC));
    if (j < 0) j = 0; if (j >= n) j = n - 1;
    den[j] += 1.;
    num[NF*j + 0] += hvar[];
    num[NF*j + 1] += zlin[];
    num[NF*j + 2] += smooth[];
    num[NF*j + 3] += curved[];
    num[NF*j + 4] += radial[];
    if (Delta < dmn[j]) dmn[j] = Delta;
    if (Delta > dmx[j]) dmx[j] = Delta;
  }

  for (int k = 0; k < NF; k++) r->err[k] = 0.;
  r->dmin = 1e30; r->dmax = -1e30; r->n = 0;
  for (int j = 0; j < n; j++) {
    if (den[j] <= 0.) continue;
    double zc = slab_z (j);
    double e[NF] = { fabs (num[NF*j+0]/den[j] - exact_hvar (zc)),
                     fabs (num[NF*j+1]/den[j] - analytic_zlin (zc)),
                     fabs (num[NF*j+2]/den[j] - exact_smooth (zc)),
                     fabs (num[NF*j+3]/den[j] - exact_curved (zc)),
                     fabs (num[NF*j+4]/den[j] - exact_radial (zc)) };
    for (int k = 0; k < NF; k++) if (e[k] > r->err[k]) r->err[k] = e[k];
    if (fabs (zc) < BAND) {
      if (dmn[j] < r->dmin) r->dmin = dmn[j];
      if (dmx[j] > r->dmax) r->dmax = dmx[j];
    }
    r->n++;
  }
}

/** Wall-clock per call, averaged over NREP. */
#ifndef NREP
# define NREP 200
#endif

double time_call (void (*fn)(Result *), Result * r)
{
  fn (r);                       // warm-up
  timer t = timer_start();
  for (int k = 0; k < NREP; k++)
    fn (r);
  return timer_elapsed (t)/NREP;
}

static void leaves_unw (Result * r) { leaves (0, r); }
static void leaves_dv  (Result * r) { leaves (1, r); }

/** Strategy 4 split: restriction alone, traversal alone. `sink` keeps them
    from being elided. */
static double sink = 0.;

static void restrict_only (Result * r)
{
  restriction (flds);
  double a = 0.;
  foreach_level (LBASE, reduction(+:a)) a += hvar[];
  r->err[0] = a;
  sink += a;
}

static void level_only (Result * r)
{
  int n = NC;
  double num[NF*n], den[n];
  for (int i = 0; i < NF*n; i++) num[i] = 0.;
  for (int i = 0; i < n; i++) den[i] = 0.;
  foreach_level (LBASE, reduction(+:num[:NF*n]) reduction(+:den[:n])) {
    int j = (int)((z - Z0)/(L0/NC));
    if (j < 0) j = 0; if (j >= n) j = n - 1;
    den[j] += 1.;
    num[NF*j + 0] += hvar[];   num[NF*j + 1] += zlin[];
    num[NF*j + 2] += smooth[]; num[NF*j + 3] += curved[];
    num[NF*j + 4] += radial[];
  }
  for (int i = 0; i < NF*n; i++) sink += num[i];
  r->n = n;
}

void report (const char * name, Result r)
{
  fprintf (stderr, "%-31s %9.2e %9.2e %9.2e %9.2e %9.2e  %8.6f %8.6f %s\n",
           name, r.err[0], r.err[1], r.err[2], r.err[3], r.err[4],
           r.dmin, r.dmax, r.dmin == r.dmax ? "uniform" : "MIXED");
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

  long nleaf = 0, nband = 0, ncyl = 0;
  foreach (reduction(+:nleaf) reduction(+:nband) reduction(+:ncyl)) {
    nleaf++;
    if (level == LBAND) nband++;
    if (level == LCYL)  ncyl++;
  }
  fprintf (stderr, "\ngrid: base %d | band %d (|z| < %g) | cylinder %d (r < %g)\n",
           LBASE, LBAND, BAND, LCYL, RCYL);
  fprintf (stderr, "leaves %ld (base-uniform %ld); at band level %ld, at cyl level %ld\n\n",
           nleaf, (long)NC*NC*NC, nband, ncyl);

  fprintf (stderr, "%-31s %9s %9s %9s %9s %9s  %8s %8s\n",
           "strategy", "hvar", "zlin", "smooth", "curved", "radial",
           "Dmin", "Dmax");
  fprintf (stderr, "%-31s %9s %9s %9s %9s %9s  %8s %8s\n",
           "", "(=1/2)", "(=z)", "(=0)", "(=tanh)", "(quad)", "", "");

  Result r;
  double t1 = time_call (leaves_unw,  &r); report ("1 foreach, unweighted", r);
  double t2 = time_call (leaves_dv,   &r); report ("2 foreach, dv()-weighted", r);
  double t3 = time_call (lattice,     &r); report ("3 foreach_region lattice", r);
  double t4 = time_call (restricted,  &r); report ("4 restriction + foreach_level", r);

  Result rr;
  double tr = time_call (restrict_only, &rr);
  double tl = time_call (level_only,    &rr);

  fprintf (stderr, "\ntiming: %d fields, %d repeats, %ld leaves\n", NF, NREP, nleaf);
  fprintf (stderr, "%-31s %10s %8s\n", "strategy", "s/call", "rel");
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "1 foreach, unweighted", t1, t1/t2);
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "2 foreach, dv()-weighted", t2, 1.0);
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "3 foreach_region lattice", t3, t3/t2);
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "4 restriction + foreach_level", t4, t4/t2);
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "  4a restriction + 1 level pass", tr, tr/t2);
  fprintf (stderr, "%-31s %10.6f %8.2f\n", "  4b foreach_level (5 fields)", tl, tl/t2);
  fprintf (stderr, "\n(sink %g)\n", sink);

  fprintf (stderr, "\nDmin/Dmax are over slabs inside the band, spanning\n"
                   "%g (base) / %g (band) / %g (cylinder)\n",
           L0/NC, L0/(1 << LBAND), L0/(1 << LCYL));
}
