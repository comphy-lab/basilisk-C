/**
# Testing that the two routes of `interfacial_power.h` agree

[interfacial_power.h](../interfacial_power.h) computes
$\Psi_\sigma = \int\mathbf{u}\cdot\mathbf{f}_\sigma\,dV$ twice over. Route 1
reads back the acceleration `a` that the solver actually applied and divides
out the face weights; route 2 rebuilds $\phi = \sigma\kappa$ and forms
$\phi\nabla f$ with `iforce.h`'s stencil. They should return the same number,
which is what makes route 2 -- the cheap one, insensitive to event ordering --
usable in place of the ground truth.

An oblate drop is released from rest with no gravity, so surface tension is
the only force in the problem and $\mathbf{a}$ contains nothing else. The
comparison is made every timestep from a `projection` event, which is the
first hook where `a` is complete and `u` is still the velocity the force acted
on: the `acceleration` events of `tension.h` and `iforce.h` have run, and
`centered.h`'s `correction()` has not. Declaring the event here rather than in
a header is what puts it first among the `projection` events, since same-named
events run in reverse declaration order.

`ag` is `{0}` because there is no gravity. Were there any, it would have to be
passed in and subtracted; and under `REDUCED` route 1 cannot be used at all,
`reduced.h` folding buoyancy into the same $\phi$.

The grid is uniform. On trees the two routes are *expected* to differ wherever
the interface sits on a refinement boundary, `iforce.h` swapping `f`'s
prolongation to the pressure's before differencing and route 2 not; that is a
separate claim and not what is tested here.

`interfacial_power_routes.asc` gets one row per timestep,
`i err_routes psi`, the deviation normalised by $|\Psi_\sigma|$ and the value
itself, so the checker can reject a run in which $\Psi_\sigma$ never departs
from zero and the agreement is vacuous. `test_interfacial_power.py` applies
the thresholds. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "acastillo/output_fields/interfacial_power.h"

#define RDROP 0.25
#define SIGMA 1.0

FILE * fp = NULL;

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;

  /**
  A density ratio, so that `alpha` varies across the interface and route 1's
  division by the face weights is actually exercised. Inviscid: viscosity
  would add nothing to $\mathbf{a}$ but would shorten the timestep. */

  rho1 = 1.;
  rho2 = 0.1;
  mu1 = mu2 = 0.;
  f.sigma = SIGMA;

#if dimension == 3
  N = 32;
#else
  N = 64;
#endif
  init_grid (N);
  run();
}

/**
An oblate spheroid: at rest it is not an equilibrium shape, so the drop
oscillates and $\Psi_\sigma$ changes sign during the run. */

event init (i = 0) {
  fraction (f, sq (RDROP) - sq (x/1.3) - sq (1.3*y) - sq (z));
  if (pid() == 0) {
    fp = fopen ("interfacial_power_routes.asc", "w");
    fprintf (fp, "# i err_routes psi\n");
  }
}

event projection (i++) {
  coord ag = {0.};
  scalar d[];
  double psi_a = interfacial_power_acceleration (ag);
  double psi_c = interfacial_power_curvature (f, d);

  /**
  Normalised by the value itself, so the tolerance means the same thing at
  every amplitude. The first timestep starts from rest and gives
  $\Psi_\sigma = 0$ exactly, which would divide by zero. */

  double scale = fabs (psi_a) > 0. ? fabs (psi_a) : 1.;
  if (pid() == 0) {
    fprintf (fp, "%d %.17g %.17g\n", i, fabs (psi_a - psi_c)/scale, psi_a);
    fflush (fp);
  }
}

event stop (i = 20) {
  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_interfacial_power_routes.ref`. */
    system ("python3 ../test_interfacial_power.py "
            "--routes interfacial_power_routes.asc 1>&2");
  }
  return 1;
}
