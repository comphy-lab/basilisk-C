/**
# `profile_*_slab` sampled at a level

`profiles_slab_restrict.h` (`restriction()` then `foreach_level()`) against
closed-form horizontal averages, on a three-level grid: base, a band about
$z=0$ at base+1, and a cylinder through the band at base+2.

Three runs of the same fields:

* **uniform** -- base level everywhere. The baseline, and the case that must
  reproduce the old leaf traversal exactly.
* **nested** -- band and cylinder refined, sampled at the base level. Must
  reproduce the uniform profile, and every slab must report one `Delta`.
* **leaves** -- the same nested grid sampled at `depth()`, which is what a
  caller reaching for the obvious default would get. Must *differ*, and must
  warn about empty slabs, or the second case proves nothing.

Fields, each isolating a different failure: `hvar` varies horizontally (slab
mean exactly 1/2) and so sees any mis-weighting between cell sizes; `zlin` is
constant within a slab and is a control that must stay exact; `smooth` has slab
mean exactly 0 and tests cancellation; `curved` is nonlinear in $z$ and
separates resampling error from weighting error; `radial` is keyed to the
cylinder, so it varies exactly where the finest cells are.

All three writers are exercised, since each restricts a different list:

* `profile_scalar_slab` -- the scalars above.
* `profile_product_slab` -- `radial*radial`, whose slab mean is
  $\overline{radial^2} = 0.1097$ against $\overline{radial}^2 = 0.0304$: a gap
  of 0.079, keyed to the refined cylinder. This is the pair that discriminates,
  since forming the product after restriction gives the product of the means.
  `zlin*smooth` rides along as a control, 0 at every height either way.

`profile_dissipation_slab` is not exercised here --
[test_dissipation_level.c](test_dissipation_level.c) covers it on a field that
is periodic in all three directions, which a velocity linear in $z$ cannot be.

Checks its own tolerances and reports through the exit code, as
[test_profile_bias.c](test_profile_bias.c) does, so it needs no `.ref`.
*/

#include "navier-stokes/centered.h"
#include "acastillo/output_fields/profiles/profiles_slab_restrict.h"

#define LBASE 5
#define LBAND (LBASE + 1)
#define LCYL  (LBASE + 2)
#define NC    (1 << LBASE)
#define BAND  (6.*L0/NC)
#define RCYL  (0.22*L0)
#define WRAMP (3.*L0/NC)
#define DTANH (3.*L0/NC)

/** The profile coordinate: y in 2D, z in 3D, as PROFILE_SLAB_COORD. */
#if dimension == 2
# define VERT y
# define HORIZ x
#else
# define VERT z
# define HORIZ x
#endif

scalar hvar[], zlin[], smooth[], curved[], radial[];

/** Distance from the vertical axis: a cylinder in 3D, a slab-centred strip in
    2D, so the innermost refinement is always off-axis and curved where it can
    be. */
double radial_profile (double px, double py) {
#if dimension == 2
  double r = fabs (px);
#else
  double r = sqrt (sq(px) + sq(py));
#endif
  return 0.5*(1. - tanh ((r - RCYL)/WRAMP));
}

/** Slab means of `radial` and `radial^2`, by quadrature over the periodic
    square. Both are height-independent: the cylinder axis is vertical. */

double radial_quad (int square)
{
  static double m1 = -1., m2 = -1.;
  if (m1 < 0.) {
    int nq = 4096;
    double h = L0/nq, s1 = 0., s2 = 0.;
    for (int i = 0; i < nq; i++)
      for (int j = 0; j < nq; j++) {
        double v = radial_profile (X0 + (i + 0.5)*h, Y0 + (j + 0.5)*h);
        s1 += v; s2 += v*v;
      }
    m1 = s1/((double) nq*nq);
    m2 = s2/((double) nq*nq);
  }
  return square ? m2 : m1;
}

/** Tolerances. `EXACT` holds the controls that are constant within a slab and
    must not move at all. `RESAMPLE` is the bound on a field the refinement
    genuinely resamples -- `curved` and `radial` -- where mode 1 differs from
    mode 0 by a discretisation error the sampling cannot remove. `DIFFER_MIN`
    is the anti-vacuity guard on mode 2. `PROD_MIN` is the gap between
    $\overline{radial^2}$ and $\overline{radial}^2$, with margin. */

#define EXACT_TOL  1e-12
#define RESAMPLE   2e-2
#define DIFFER_MIN 5e-2
#define PROD_MIN   5e-2

int nfail = 0;
double got_prod[3];   // slab mean of radial*radial, per mode

int mode;   // 0 uniform, 1 nested@base, 2 nested@leaves

/** Column pairs of the scalar block, one (mean, mean^2) per field. */
#define COL_HVAR   3
#define COL_ZLIN   5
#define COL_SMOOTH 7
#define COL_CURVED 9
#define COL_RADIAL 11
#define NCOL       13

/** Max |mode a - mode b| over slabs, in one column of the scalar block. */

double column_gap (int ma, int mb, int col)
{
  char na[80], nb[80];
  sprintf (na, "prof_level_mode%d.asc", ma);
  sprintf (nb, "prof_level_mode%d.asc", mb);
  FILE * fa = fopen (na, "r"), * fb = fopen (nb, "r");
  if (fa == NULL || fb == NULL) { perror (na); exit (1); }
  char la[2048], lb[2048];
  while (fgets (la, sizeof(la), fa) && la[0] == '#');
  while (fgets (lb, sizeof(lb), fb) && lb[0] == '#');
  double gap = 0.;
  do {
    double va[NCOL], vb[NCOL];
    int ca = 0, cb = 0;
    char * pa = la, * pb = lb; int nch;
    while (ca < NCOL && sscanf (pa, "%lf%n", &va[ca], &nch) == 1) { pa += nch; ca++; }
    while (cb < NCOL && sscanf (pb, "%lf%n", &vb[cb], &nch) == 1) { pb += nch; cb++; }
    if (ca > col && cb > col) {
      double d = fabs (va[col] - vb[col]);
      if (d > gap) gap = d;
    }
  } while (fgets (la, sizeof(la), fa) && fgets (lb, sizeof(lb), fb));
  fclose (fa); fclose (fb);
  return gap;
}

int main() {
  L0 = 1.; X0 = Y0 = Z0 = -L0/2.;
  for (mode = 0; mode <= 2; mode++) {
    init_grid (NC);
    run();
  }

  if (pid() != 0)
    return 0;

  /** The one true control: `zlin` is constant within a slab, so every cell of
      a slab holds the same value whatever its size and the sampling level
      cannot touch it. `smooth` varies horizontally and is resampled like
      `curved` and `radial`, so it belongs with those. */
  double gz = column_gap (0, 1, COL_ZLIN);
  fprintf (stderr, "control zlin    nested@base vs uniform : %.3e\n", gz);
  if (gz > EXACT_TOL) {
    fprintf (stderr, "FAIL: control zlin moved; the binning itself is wrong\n");
    nfail++;
  }

  /** `hvar` is exactly 1/2 in every slab whatever the grid. */
  double gh = column_gap (0, 1, COL_HVAR);
  fprintf (stderr, "hvar    nested@base vs uniform : %.3e\n", gh);
  if (gh > EXACT_TOL) {
    fprintf (stderr, "FAIL: hvar moved across the level jump\n");
    nfail++;
  }

  /** Resampled fields: mode 1 keeps a discretisation residual, bounded. */
  struct { int col; const char * name; } res[] = {
    {COL_SMOOTH, "smooth"}, {COL_CURVED, "curved"}, {COL_RADIAL, "radial"}
  };
  for (int k = 0; k < 3; k++) {
    double g = column_gap (0, 1, res[k].col);
    fprintf (stderr, "%-7s nested@base vs uniform : %.3e  (tol %g)\n",
             res[k].name, g, RESAMPLE);
    if (g > RESAMPLE) {
      fprintf (stderr, "FAIL: %s residual is larger than resampling alone\n",
               res[k].name);
      nfail++;
    }
  }

  /** Anti-vacuity: sampling the same grid at the leaves must visibly differ. */
  double gv = column_gap (0, 2, COL_RADIAL);
  fprintf (stderr, "radial  leaves vs uniform      : %.3e  (min %g)\n",
           gv, DIFFER_MIN);
  if (gv < DIFFER_MIN) {
    fprintf (stderr, "FAIL: sampling at the leaves does not differ: "
             "the test is vacuous\n");
    nfail++;
  }

  /** The product is formed at the leaves, so it is the mean of the square,
      not the square of the mean. */
  double ex2 = radial_quad (1), ex1 = radial_quad (0);
  fprintf (stderr, "\nmean(radial*radial) = %.6f  (exact %.6f, "
           "square-of-mean %.6f)\n", got_prod[1], ex2, ex1*ex1);
  if (fabs (got_prod[1] - ex2) > RESAMPLE) {
    fprintf (stderr, "FAIL: product is not the mean of the square\n");
    nfail++;
  }
  if (fabs (ex2 - ex1*ex1) < PROD_MIN) {
    fprintf (stderr, "FAIL: the two forms are too close to tell apart\n");
    nfail++;
  }

  if (nfail) {
    fprintf (stderr, "FAIL: %d check(s) failed.\n", nfail);
    return 1;
  }
  fprintf (stderr, "PASS: profiles sample correctly at slablevel.\n");
  return 0;
}

event init (i = 0) {
  /** `refine()` is a tree operation: on multigrid the hierarchy is complete and
      uniform, so there is nothing to nest and modes 1/2 degenerate to mode 0.
      Kept compiling there because the profile routines themselves are grid
      agnostic -- restriction and foreach_level exist on multigrid too. */
#if TREE
  if (mode > 0) {
    refine (fabs(VERT) < BAND && level < LBAND);
# if dimension == 2
    refine (fabs(VERT) < BAND && sq(x) < sq(RCYL) && level < LCYL);
# else
    refine (fabs(VERT) < BAND && sq(x) + sq(y) < sq(RCYL) && level < LCYL);
# endif
  }
#endif
  foreach() {
    hvar[]   = (HORIZ < 0.) ? 1. : 0.;
    zlin[]   = VERT;
#if dimension == 2
    smooth[] = cos (2.*pi*x/L0);
#else
    smooth[] = cos (2.*pi*x/L0)*cos (2.*pi*y/L0);
#endif
    curved[] = tanh (VERT/DTANH);
    radial[] = radial_profile (x, y);
    /** Linear in space: every velocity gradient is an exact constant, so the
        dissipation profile has a closed form even across a level jump. */
    u.x[] = 2.*VERT; // only the vertical may vary linearly: x,y are periodic
    u.y[] = 0.;
#if dimension == 3
    u.z[] = 0.;
#endif
  }
  boundary ((scalar *){u});
}

event go (i = 0) {
  char name[80];
  sprintf (name, "prof_level_mode%d.asc", mode);
  double del = L0/NC;
  /** On a decomposed multigrid each rank owns a sub-box, so `depth()` is
      rank-local and smaller than the global LBASE: asking for a level the local
      tree does not have reads past its arrays. Clamp. */
  int slev = (mode == 2) ? depth() : min (LBASE, depth());
  double vmin = (dimension == 2 ? Y0 : Z0);

  scalar * list = {hvar, zlin, smooth, curved, radial};
  profile_scalar_slab (list, filename = name,
                       hmin = vmin + del/2., hmax = vmin + L0 - del/2.,
                       n = NC, slablevel = slev, mode = "w");

  /** `radial*radial` discriminates: forming the product after restriction
      would give the product of the means, 0.0304 rather than 0.1097.
      `zlin*smooth` is a control, 0 at every height either way. */
  char pname[80]; sprintf (pname, "prod_level_mode%d.asc", mode);
  scalar * l1 = {radial, zlin};
  scalar * l2 = {radial, smooth};
  profile_product_slab (l1, l2, filename = pname,
                        hmin = vmin + del/2., hmax = vmin + L0 - del/2.,
                        n = NC, slablevel = slev, mode = "w");

  if (pid() == 0) {
    fprintf (stderr, "mode %d: sampled at level %d, depth %d\n",
             mode, slev, depth());

    /** Read the product block back: column 3 is mean(radial*radial). */
    FILE * fp = fopen (pname, "r");
    if (fp == NULL) { perror (pname); exit (1); }
    char line[1024];
    while (fgets (line, sizeof(line), fp) && line[0] == '#');
    double sum = 0.; int nrow = 0, iprof;
    double z, delta, p1, p2;
    do {
      if (sscanf (line, "%d %lf %lf %lf %lf",
                  &iprof, &z, &delta, &p1, &p2) == 5) {
        sum += p1; nrow++;
      }
    } while (fgets (line, sizeof(line), fp));
    fclose (fp);
    got_prod[mode] = nrow ? sum/nrow : 0.;
    fprintf (stderr, "        mean(radial*radial) = %.6f\n", got_prod[mode]);
  }
}

event stop (i = 1) { return 1; }
