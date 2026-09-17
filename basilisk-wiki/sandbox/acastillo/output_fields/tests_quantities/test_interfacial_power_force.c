/**
# `surface_tension_force()` reproduces `interfacial_power_curvature()`'s decomposition

`interfacial_power_curvature()` is now built on `surface_tension_force()`: the
vector $\mathbf{F}_\sigma = \sigma\kappa\nabla f$, from which
`u`.$\mathbf{F}_\sigma$ integrates to $\Psi_\sigma$. This checks that
decomposition end to end -- pointwise and integrated -- on an off-centre drop
with an arbitrary (non-uniform, non-radial) velocity field, so the check is
not accidentally satisfied by one of the symmetric cases
[test_interfacial_power_static.c](test_interfacial_power_static.c) already
covers.

Both are exact identities: the same $\phi$ field and the same stencil,
computed once inside `surface_tension_force()` and reused by
`interfacial_power_curvature()`. So `ROUNDOFF_TOL` applies directly, with no
order-of-accuracy argument -- this isn't testing `curvature()`'s accuracy,
only that the new primitive is a faithful decomposition of the existing one.
`fraction()` is accurate enough for that; `Vofi` is not needed here. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "acastillo/output_fields/interfacial_power.h"

#define RDROP 0.25
#define SIGMA 1.3

int main() {
  L0 = 1.;
  X0 = Y0 = Z0 = -L0/2.;
  f.sigma = SIGMA;

#if dimension == 3
  N = 32;
#else
  N = 64;
#endif
  init_grid (N);

  // Off-centre, so the check is not vacuously satisfied by mirror symmetry.
#if dimension == 3
  fraction (f, sq (RDROP) - sq (x - 0.07) - sq (y + 0.043) - sq (z - 0.031));
#else
  fraction (f, sq (RDROP) - sq (x - 0.07) - sq (y + 0.043));
#endif

  // Neither uniform nor radial, so this exercises a case the closed-form
  // checks in test_interfacial_power_static.c cannot.
  foreach() {
    u.x[] = sin (5.*x) + cos (3.*y);
    u.y[] = cos (4.*x) - sin (2.*y);
#if dimension == 3
    u.x[] += cos (2.*z);
    u.y[] += sin (3.*z);
    u.z[] = sin (x + y + z);
#endif
  }

  vector Fs[];
  surface_tension_force (f, Fs);

  scalar d[];
  double psi_c = interfacial_power_curvature (f, d);

  double dmax = 0., err_field = 0., psi_dot = 0.;
  foreach (reduction(max:dmax) reduction(max:err_field)
          reduction(+:psi_dot)) {
    double di = 0.;
    foreach_dimension()
      di += u.x[]*Fs.x[];
    dmax = max (dmax, fabs (d[]));
    err_field = max (err_field, fabs (di - d[]));
    psi_dot += dv()*di;
  }

  double scale = fabs (psi_c) > 0. ? fabs (psi_c) : 1.;
  double fscale = dmax > 0. ? dmax : 1.;

  if (pid() == 0) {
    FILE * fp = fopen ("interfacial_power_force.asc", "w");
    fprintf (fp, "# N err_field err_psi psi\n");
    fprintf (fp, "%d %.17g %.17g %.17g\n", N, err_field/fscale,
             fabs (psi_dot - psi_c)/scale, psi_c);
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_interfacial_power_force.ref`. */
    system ("python3 ../test_interfacial_power.py "
           "--force interfacial_power_force.asc 1>&2");
  }
}
