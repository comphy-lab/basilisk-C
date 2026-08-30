/**
# Testing the `profiles/` routines on a two-phase shear flow

In this example use the sample results to test `profile_foreach_region()`,
`profile_product_foreach_region()` and `profile_dissipation_foreach_region()`
on the same field: a two-phase shear flow, split by a horizontal interface
at y = 0:

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
  profile_foreach_region(list, unity, "profiles_scalar.asc", SAMPLING);

  scalar * list1 = {u.x, u.y};
  scalar * list2 = {u.y, u.y};
  profile_product_foreach_region(list1, list2, unity, "profiles_product.asc", SAMPLING);

  profile_dissipation_foreach_region(u, unity, "profiles_dissipation.asc", SAMPLING);

  #undef SAMPLING

  /**
  Verify the three profiles that were just written: the means/variances,
  cross products, and velocity-gradient/strain/dissipation terms must all
  match the closed form of this affine two-phase field on every plane. The
  report goes to stderr, i.e. to the `log` diffed against
  `test_profiles.ref`. */
  if (pid() == 0)
    system ("python3 ../test_profiles.py --scalar profiles_scalar.asc "
            "--product profiles_product.asc "
            "--dissipation profiles_dissipation.asc 1>&2");
}
