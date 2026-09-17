/**
# `profiles_slab.h` against `profiles.h`, where they must agree

[profiles_slab.h](../profiles/profiles_slab.h) bins each cell into a slab in
one traversal where [profiles.h](../profiles/profiles.h) point-samples every
slab in turn. The header claims the two coincide when the sampled region is at
one resolution and each slab is exactly one row of cells. This runs both on
the same field under exactly that geometry and compares the files they write,
column by column.

The field is the two-phase shear flow of [test_profiles.c](test_profiles.c):

  f = 1 for y < 0 (phase 1), f = 0 for y > 0 (phase 2)
  rho = rho1*f + rho2*(1-f), mu = mu1*f + mu2*(1-f)

  U_i(x,y) = alpha_shear*x + beta_i*y
  V_i(x,y) = -alpha_shear*y

with mu1*beta1 = mu2*beta2. There is no refinement here: the grid is uniform
even on `quadtree`/`octree`, which is the point -- the precondition is about
resolution, not about the grid type.

The sampling geometry is what makes the claim testable. With `del = L0/N`, the
slab centres are the cell centres, `hmin = Y0 + del/2`, `hmax = Y0 + L0 -
del/2` and `n = N`, so `deltahn` comes out to `del` and slab `j` holds cell row
`j` alone. The in-plane sample counts `m1` (and `m2`) are `N` over the same
inset range, so `foreach_region()` lands one sample on each cell centre.
Getting this wrong is the interesting failure: see
[test_profiles_slab_bias.c](test_profiles_slab_bias.c).

The two accumulate in different orders -- one over sampled planes, the other
over cells in grid order -- so `test_profiles_slab.py` compares them against a
tolerance rather than asserting bit-for-bit equality. In practice they come out
identical here, bar a column of the octree dissipation profile that is
analytically zero and so holds only roundoff. */

#define LEVEL 5

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "acastillo/output_fields/profiles/profiles.h"
#include "acastillo/output_fields/profiles/profiles_slab.h"

double tau = 1.;
double alpha_shear = 0.5;

int main(){

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid (N);

  rho1 = 3.; rho2 = 1.;
  mu1 = 3.; mu2 = 1.;
  double beta1 = tau / mu1;
  double beta2 = tau / mu2;

#if dimension == 3
  fraction (f, -z);
#else
  fraction (f, -y);
#endif

  foreach() {
    double beta = f[] ? beta1 : beta2;
    u.x[] = alpha_shear * x + beta * y;
    u.y[] = -alpha_shear * y;
#if dimension == 3
    u.z[] = alpha_shear * z + beta * x;
#endif
  }
  boundary ((scalar *){u});
  foreach()
    rhov[] = rho1 * f[] + rho2 * (1. - f[]);
  boundary ({rhov});

  /**
  One slab per cell row, one in-plane sample per cell centre. */
  double del = L0/N;
  #define SAMPLING \
    xmin = X0 + del/2., xmax = X0 + L0 - del/2., \
    ymin = Y0 + del/2., ymax = Y0 + L0 - del/2., \
    hmin = Z0 + del/2., hmax = Z0 + L0 - del/2., \
    n = N, m1 = N, m2 = N, mode = "w"

  scalar * list = {u.x, u.y, rhov};
  profile_foreach_region (list, unity, "region_scalar.asc", SAMPLING);
  profile_scalar_slab    (list, unity, "slab_scalar.asc",   SAMPLING);

  scalar * list1 = {u.x, u.y};
  scalar * list2 = {u.y, u.y};
  profile_product_foreach_region (list1, list2, unity, "region_product.asc",
                                  SAMPLING);
  profile_product_slab           (list1, list2, unity, "slab_product.asc",
                                  SAMPLING);

  profile_dissipation_foreach_region (u, unity, "region_dissipation.asc",
                                      SAMPLING);
  profile_dissipation_slab           (u, unity, "slab_dissipation.asc",
                                      SAMPLING);

  #undef SAMPLING

  /**
  The report goes to stderr, i.e. to the `log` diffed against
  `test_profiles_slab.ref`. */
  if (pid() == 0)
    system ("python3 ../test_profiles_slab.py "
            "--pair region_scalar.asc slab_scalar.asc "
            "--pair region_product.asc slab_product.asc "
            "--pair region_dissipation.asc slab_dissipation.asc 1>&2");
}
