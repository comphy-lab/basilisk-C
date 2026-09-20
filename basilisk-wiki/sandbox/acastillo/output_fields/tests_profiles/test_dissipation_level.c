/**
# Dissipation profiles on a mixed-resolution grid

`profile_dissipation_slab` against closed-form slab means, on the three grids of
[test_profiles_level.c](test_profiles_level.c): uniform, nested sampled at the
base level, nested sampled at the leaves.

Two velocity fields, both periodic in all three directions so every stencil is
exact and no cell has to be excluded -- unlike
[test_strain_and_vorticity_smooth.c](../../../../wiki/sandbox/acastillo/output_fields/tests_quantities/test_strain_and_vorticity_smooth.c),
which insets one layer because it declares no periodicity. A slab average cannot
inset: the end slabs would have no cells left.

* **abc** -- ABC flow, $k = 2\pi/L_0$. Beltrami, so $\omega = k u$ and
  $\Omega^2 = k^2|u|^2/2$; the horizontal average of $|u|^2$ is $A^2+B^2+C^2$ at
  every height, giving $\overline{S^2} = \overline{\Omega^2} =
  k^2(A^2+B^2+C^2)/2$, constant in $z$. The reference comes from the Beltrami
  identity, not from a transcription of the code's expression.
* **strain** -- $u = \nabla\phi$ with $\phi = \cos(kx)\cos(ky)/k$, irrotational
  by construction: $\overline{S^2} = k^2$ and $\overline{\Omega^2} = 0$
  **exactly**. ABC alone cannot catch a bug that swaps $S$ for $\Omega$, since
  it has $S^2 = \Omega^2$; this field fails loudly if one leaks into the other.

Both are $z$-independent in the mean, so every slab has the same answer and a
level jump shows up as a departure from a constant rather than as a shape
change. `dxx` and `dzz` average to 0 in both fields.

Checks its own tolerances and reports through the exit code, as
[test_profile_bias.c](test_profile_bias.c) does, so it needs no `.ref`. The
tolerances scale with $\Delta^2$: at $k = 2\pi$, $N = 32$ the second-order
error is ~1.3%, not a roundoff bound.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "acastillo/output_fields/profiles/profiles_slab_restrict.h"

#define LBASE 5
#define LBAND (LBASE + 1)
#define LCYL  (LBASE + 2)
#define NC    (1 << LBASE)
#define BAND  (6.*L0/NC)
#define RCYL  (0.22*L0)

double kw, Aabc = 1., Babc = 2., Cabc = 3.;

int mode;    // 0 uniform, 1 nested@base, 2 nested@leaves
int field;   // 0 abc, 1 strain

/** Tolerances. `DISS_REL` is the O(Delta^2) bound the closed forms are met to;
    `JUMP` is the result being tested -- sampling the nested grid at the base
    level must cost nothing; `DIFFER` is the anti-vacuity guard, since mode 2
    samples the same grid at the leaves and must visibly bias it. Omega^2 on
    the irrotational field is zero to roundoff on a uniform grid but not across
    a level jump: prolongation into the refined band is not curl-free, and the
    two slabs straddling each band edge carry a spurious Omega^2 (measured
    6.3e-3, against 1.6e-5 inside the band and ~1e-30 outside). */

#define DISS_REL   0.05
#define JUMP_TOL   0.01
#define DIFFER_MIN 0.2
#define IRROT_TOL  1e-12
#define IRROT_JUMP 2e-2
#define FLAT_TOL   0.05

/** Slab means of S2 and O2 per mode, filled by `go` and checked in `main`. */
double got_S2[3], got_O2[3], spread_S2[3], max_O2[3], mean_dzz[3];
int nfail = 0;

/** Closed-form slab means, both constant in z. */

double exact_S2 (void)
{
  return field == 0 ? 0.5*sq(kw)*(sq(Aabc) + sq(Babc) + sq(Cabc)) : sq(kw);
}

double exact_O2 (void)
{
  return field == 0 ? 0.5*sq(kw)*(sq(Aabc) + sq(Babc) + sq(Cabc)) : 0.;
}

/** All three modes of one field are in; check them. */

static void check_field (void)
{
  double ex_s2 = exact_S2(), ex_o2 = exact_O2();
  double scale = max (ex_s2, 1.);

  /** Modes 0 and 1 must both land on the closed form. */
  for (int m = 0; m <= 1; m++) {
    double e = fabs (got_S2[m] - ex_s2)/scale;
    if (e > DISS_REL) {
      fprintf (stderr, "FAIL: field %d mode %d S2 %.4f vs exact %.4f "
               "(rel %.3e > %g)\n", field, m, got_S2[m], ex_s2, e, DISS_REL);
      nfail++;
    }
  }

  /** The result: the level jump must not move the profile. */
  double d = fabs (got_S2[0] - got_S2[1])/scale;
  if (d > JUMP_TOL) {
    fprintf (stderr, "FAIL: field %d nested@base differs from uniform "
             "(rel %.3e > %g)\n", field, d, JUMP_TOL);
    nfail++;
  }

  /** Anti-vacuity: sampling the same grid at the leaves must bias it. */
  double d2 = fabs (got_S2[0] - got_S2[2])/scale;
  if (d2 < DIFFER_MIN) {
    fprintf (stderr, "FAIL: field %d sampling at the leaves does not differ "
             "(rel %.3e < %g): the test is vacuous\n", field, d2, DIFFER_MIN);
    nfail++;
  }

  /** Omega^2 separates the symmetric and antisymmetric parts. */
  if (ex_o2 == 0.) {
    if (max_O2[0] > IRROT_TOL) {
      fprintf (stderr, "FAIL: field %d Omega^2 %.3e on a uniform grid is not "
               "zero (> %g)\n", field, max_O2[0], IRROT_TOL);
      nfail++;
    }
    if (max_O2[1] > IRROT_JUMP) {
      fprintf (stderr, "FAIL: field %d spurious Omega^2 %.3e at the level jump "
               "(> %g)\n", field, max_O2[1], IRROT_JUMP);
      nfail++;
    }
  }
  else {
    double e = fabs (got_O2[1] - ex_o2)/ex_o2;
    if (e > DISS_REL) {
      fprintf (stderr, "FAIL: field %d Omega^2 %.4f vs exact %.4f "
               "(rel %.3e > %g)\n", field, got_O2[1], ex_o2, e, DISS_REL);
      nfail++;
    }
  }

  /** The exact means are z-independent, so a jump shows up as a bump. */
  for (int m = 0; m <= 1; m++)
    if (spread_S2[m]/scale > FLAT_TOL) {
      fprintf (stderr, "FAIL: field %d mode %d S2 varies with z "
               "(spread %.3e > %g)\n", field, m, spread_S2[m]/scale, FLAT_TOL);
      nfail++;
    }

  /** dzz averages to zero in both fields. */
  for (int m = 0; m <= 1; m++)
    if (fabs (mean_dzz[m]) > IRROT_TOL) {
      fprintf (stderr, "FAIL: field %d mode %d dzz %.3e does not average to "
               "zero\n", field, m, mean_dzz[m]);
      nfail++;
    }
}

int main()
{
  L0 = 1.; X0 = Y0 = Z0 = -L0/2.;
  kw = 2.*pi/L0;
  periodic (right); periodic (top); periodic (front);
  for (field = 0; field <= 1; field++) {
    for (mode = 0; mode <= 2; mode++) {
      init_grid (NC);
      run();
    }
    if (pid() == 0)
      check_field();
  }

  if (pid() == 0) {
    if (nfail) {
      fprintf (stderr, "FAIL: %d check(s) failed.\n", nfail);
      return 1;
    }
    fprintf (stderr, "PASS: dissipation profiles sample correctly at "
             "slablevel across a level jump.\n");
  }
  return 0;
}

event init (i = 0)
{
#if TREE
  if (mode > 0) {
    refine (fabs(z) < BAND && level < LBAND);
    refine (fabs(z) < BAND && sq(x) + sq(y) < sq(RCYL) && level < LCYL);
  }
#endif
  foreach() {
    if (field == 0) {
      u.x[] = Aabc*sin(kw*z) + Cabc*cos(kw*y);
      u.y[] = Babc*sin(kw*x) + Aabc*cos(kw*z);
      u.z[] = Cabc*sin(kw*y) + Babc*cos(kw*x);
    }
    else {
      u.x[] = -sin(kw*x)*cos(kw*y);
      u.y[] = -cos(kw*x)*sin(kw*y);
      u.z[] = 0.;
    }
  }
  boundary ((scalar *){u});
}

event go (i = 0)
{
  char name[80];
  sprintf (name, "dissipation_level_mode%d_field%d.asc", mode, field);
  double del = L0/NC;
  /** `depth()` is rank-local under MPI, so clamp as test_profiles_level does. */
  int slev = (mode == 2) ? depth() : min (LBASE, depth());

  profile_dissipation_slab (u, filename = name,
                            hmin = Z0 + del/2., hmax = Z0 + L0 - del/2.,
                            n = NC, slablevel = slev, mode = "w");

  /** Read the block back and reduce it, as test_profile_bias.c does. */
  if (pid() == 0) {
    FILE * fp = fopen (name, "r");
    if (fp == NULL) { perror (name); exit (1); }
    char line[1024];
    while (fgets (line, sizeof(line), fp) && line[0] == '#');

    double s2sum = 0., o2sum = 0., dzzsum = 0.;
    double s2min = HUGE_VAL, s2max = -HUGE_VAL, o2absmax = 0.;
    int nrow = 0, iprof;
    double z, delta, s2, s2sq, o2, o2sq, nus2, nus2sq, dxx, dxxsq, dzz, dzzsq;
    do {
      if (sscanf (line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &iprof, &z, &delta, &s2, &s2sq, &o2, &o2sq,
                  &nus2, &nus2sq, &dxx, &dxxsq, &dzz, &dzzsq) == 13) {
        s2sum += s2; o2sum += o2; dzzsum += dzz;
        if (s2 < s2min) s2min = s2;
        if (s2 > s2max) s2max = s2;
        if (fabs(o2) > o2absmax) o2absmax = fabs(o2);
        nrow++;
      }
    } while (fgets (line, sizeof(line), fp));
    fclose (fp);

    if (nrow == 0) {
      fprintf (stderr, "FAIL: no data rows in %s\n", name);
      nfail++;
      return 0;
    }
    got_S2[mode]    = s2sum/nrow;
    got_O2[mode]    = o2sum/nrow;
    mean_dzz[mode]  = dzzsum/nrow;
    max_O2[mode]    = o2absmax;
    spread_S2[mode] = s2max - s2min;

    fprintf (stderr, "field %d mode %d (level %d): S2 %10.4f  O2 %10.4f\n",
             field, mode, slev, got_S2[mode], got_O2[mode]);
  }
}

