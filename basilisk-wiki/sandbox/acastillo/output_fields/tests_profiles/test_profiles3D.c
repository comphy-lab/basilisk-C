/**
# Testing the `profiles/` routines on a two-phase shear flow, in 3D

3D analogue of `test_profiles.c`. In 3D, `profile_*_foreach_region()` sweeps
the x,y plane and profiles along z (`SETUP_PROFILE_PLANE` samples `hprof` as
the third coordinate), so the field and interface are relabeled accordingly:
a two-phase shear flow split by a horizontal interface at z = 0.

  f = 1 for z < 0 (phase 1), f = 0 for z > 0 (phase 2)
  rho = rho1*f + rho2*(1-f), mu = mu1*f + mu2*(1-f)

  U_i(x,y,z) = alpha_shear*x + beta_i*z
  V_i(x,y,z) = 0
  W_i(x,y,z) = -alpha_shear*z

with beta_i chosen so the shear stress mu_i*beta_i is continuous across the
interface (tau = mu1*beta1 = mu2*beta2), even though beta_i itself jumps.
As in 2D, div(u) = dU/dx + dV/dy + dW/dz = alpha_shear + 0 - alpha_shear = 0.

*/

#define LEVEL 5
#define MAXLEVEL 6

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "acastillo/output_fields/profiles/profiles.h"

double tau = 1.;
double alpha_shear = 0.5;

int main(){

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  #if TREE
    // Refine near the interface at z = 0, where the field is discontinuous.
    double band = 4. * L0 / (1 << LEVEL);
    refine (fabs(z) < band && level < MAXLEVEL);
  #endif

  rho1 = 3.; rho2 = 1.;
  mu1 = 3.; mu2 = 1.;
  double beta1 = tau / mu1;
  double beta2 = tau / mu2;

  fraction (f, -z);

  foreach() {
    double beta = f[] ? beta1 : beta2;
    u.x[] = alpha_shear * x + beta * z;
    u.y[] = 0.;
    u.z[] = -alpha_shear * z;
  }
  boundary({u.x, u.y, u.z});
  foreach()
    rhov[] = rho1 * f[] + rho2 * (1. - f[]);
  boundary({rhov});

  // Sampling geometry shared by every profile call below: same domain,
  // same plane resolution, so every routine samples the field identically.
  // The swept plane is x,y (inset from the domain edges so the dissipation
  // stencil's neighbor access stays inside the domain); the profile
  // coordinate is z (also inset, same reasoning as the 2D test).
  #define SAMPLING \
    xmin = X0 + L0 / N, xmax = X0 + L0 - L0 / N, \
    ymin = Y0 + L0 / N, ymax = Y0 + L0 - L0 / N, \
    hmin = Z0 + L0 / N, hmax = Z0 + L0 - L0 / N, \
    n = 16, m1 = 32, m2 = 32, mode = "w"

  scalar * list = {u.x, u.y, u.z, rhov};
  profile_foreach_region(list, unity, "profiles_scalar3D.asc", SAMPLING);

  scalar * list1 = {u.x, u.z};
  scalar * list2 = {u.z, u.z};
  profile_product_foreach_region(list1, list2, unity, "profiles_product3D.asc", SAMPLING);

  profile_dissipation_foreach_region(u, unity, "profiles_dissipation3D.asc", SAMPLING);

  #undef SAMPLING

  /**
  Verify the three profiles that were just written: the means/variances,
  cross products, and velocity-gradient/strain/dissipation terms must all
  match the closed form of this affine two-phase field on every plane. The
  report goes to stderr, i.e. to the `log` diffed against
  `test_profiles3D.ref`. */
  if (pid() == 0)
    system ("python3 ../test_profiles.py --scalar3D profiles_scalar3D.asc "
            "--product3D profiles_product3D.asc "
            "--dissipation3D profiles_dissipation3D.asc 1>&2");
}
