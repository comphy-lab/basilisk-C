/**
# Testing the `profiles/` routines on a two-phase shear flow, weighted by rho

Same field and setup as `test_profiles.c`, but every profile is weighted by
`rhov` instead of `unity`. Since the interface is horizontal, `rhov` is
constant on any single sampled plane (rho1 for y<0, rho2 for y>0), so the
weighted mean equals the unweighted one exactly:

  <Q>_w(y) = sum(rho(y)*Q(x,y)) / sum(rho(y)) = sum(Q(x,y)) / n = <Q>(y)

This does not change the expected numbers -- the closed forms are identical
to `test_profiles.c`'s -- but it does exercise the weighted code path
(`weight = w[]`, normalized by `total_weight`) in every routine, rather than
the `w.i == unity.i` shortcut that skips it.

  f = 1 for y < 0 (phase 1), f = 0 for y > 0 (phase 2)
  rho = rho1*f + rho2*(1-f), mu = mu1*f + mu2*(1-f)

  U_i(x,y) = alpha_shear*x + beta_i*y
  V_i(x,y) = -alpha_shear*y

with beta_i chosen so the shear stress mu_i*beta_i is continuous across the
interface (tau = mu1*beta1 = mu2*beta2), even though beta_i itself jumps.

*/

#define LEVEL 6
#define MAXLEVEL 8

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "acastillo/output_fields/profiles/profiles.h"

double tau = 1.;
double alpha_shear = 0.5;

int main(){

  L0 = 1.0;
  X0 = Y0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  #if TREE
    // Refine near the interface at y = 0, where the field is discontinuous.
    double band = 4. * L0 / (1 << LEVEL);
    refine (fabs(y) < band && level < MAXLEVEL);
  #endif

  rho1 = 3.; rho2 = 1.;
  mu1 = 3.; mu2 = 1.;
  double beta1 = tau / mu1;
  double beta2 = tau / mu2;

  fraction (f, -y);

  foreach() {
    double beta = f[] ? beta1 : beta2;
    u.x[] = alpha_shear * x + beta * y;
    u.y[] = -alpha_shear * y;
  }
  boundary({u.x, u.y});
  foreach()
    rhov[] = rho1 * f[] + rho2 * (1. - f[]);
  boundary({rhov});

  // Sampling geometry shared by every profile call below: same domain,
  // same plane resolution, so every routine samples the field identically.
  #define SAMPLING \
    xmin = X0 + L0 / N, xmax = X0 + L0 - L0 / N, \
    hmin = Y0 + L0 / N, hmax = Y0 + L0 - L0 / N, \
    n = 32, m1 = 128, mode = "w"

  scalar * list = {u.x, u.y, rhov};
  profile_foreach_region(list, rhov, "profiles_scalar_weighted.asc", SAMPLING);

  scalar * list1 = {u.x, u.y};
  scalar * list2 = {u.y, u.y};
  profile_product_foreach_region(list1, list2, rhov, "profiles_product_weighted.asc", SAMPLING);

  profile_dissipation_foreach_region(u, rhov, "profiles_dissipation_weighted.asc", SAMPLING);

  #undef SAMPLING

  /**
  Verify the three profiles that were just written: the rho-weighted
  means/variances, cross products, and velocity-gradient/strain/dissipation
  terms must all match the same closed form as the unweighted case, since
  rho is constant on every horizontal plane. The report goes to stderr, i.e.
  to the `log` diffed against `test_profiles_weighted.ref`. */
  if (pid() == 0)
    system ("python3 ../test_profiles.py --scalar profiles_scalar_weighted.asc "
            "--product profiles_product_weighted.asc "
            "--dissipation profiles_dissipation_weighted.asc 1>&2");
}
