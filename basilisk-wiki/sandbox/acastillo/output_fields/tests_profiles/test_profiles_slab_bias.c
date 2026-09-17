/**
# `profiles_slab.h` against `profiles.h`, where they must not agree

The companion to [test_profiles_slab.c](test_profiles_slab.c), which pins the
agreement under the documented geometry. This pins the other half of the
claim: that the precondition is load-bearing rather than defensive. The field,
the sampling geometry and the calls are the same; the one change is a band
refined one level around the interface, so the sampled region no longer sits
at a single resolution.

There the two coincide because `foreach_region()` puts one sample on each cell
centre. Here a slab that crosses the refined band contains two rows of fine
cells, and the two routines part company:
[profiles.h](../profiles/profiles.h) still takes one plane of samples through
the slab, while [profiles_slab.h](../profiles/profiles_slab.h) bins every cell
in it. Neither is wrong -- they are different averages -- but they are not
interchangeable, which is what the header says and what this measures.

Tree grids only: on `multigrid` the `refine()` below is skipped, the grid stays
uniform, and the test would be [test_profiles_slab.c](test_profiles_slab.c)
again under another name.

`test_profiles_slab.py` checks the deviation is large, and fails if the two
ever start agreeing -- which would mean `profiles_slab.h` had outgrown its
documentation. */

#define LEVEL 5
#define MAXLEVEL 6

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

  /**
  The one difference from `test_profiles_slab.c`: two resolutions inside the
  sampled region. */
#if TREE
  double band = 4. * L0 / (1 << LEVEL);
  #if dimension == 3
    refine (fabs(z) < band && level < MAXLEVEL);
  #else
    refine (fabs(y) < band && level < MAXLEVEL);
  #endif
#endif

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

  double del = L0/N;
  #define SAMPLING \
    xmin = X0 + del/2., xmax = X0 + L0 - del/2., \
    ymin = Y0 + del/2., ymax = Y0 + L0 - del/2., \
    hmin = Z0 + del/2., hmax = Z0 + L0 - del/2., \
    n = N, m1 = N, m2 = N, mode = "w"

  scalar * list = {u.x, u.y, rhov};
  profile_foreach_region (list, unity, "region_scalar.asc", SAMPLING);
  profile_scalar_slab    (list, unity, "slab_scalar.asc",   SAMPLING);

  #undef SAMPLING

  /**
  The report goes to stderr, i.e. to the `log` diffed against
  `test_profiles_slab_bias.ref`. */
  if (pid() == 0)
    system ("python3 ../test_profiles_slab.py "
            "--differ region_scalar.asc slab_scalar.asc 1>&2");
}
