/**
# Testing $\Psi_\sigma = -d(\sigma A)/dt$ on a relaxing drop

[interfacial_power.h](../interfacial_power.h) rests on the identity

$$
\Psi_\sigma = -\frac{d(\sigma A)}{dt}
$$

exact in the continuum, since $\phi\,\mathbf{u}\cdot\nabla f =
\nabla\cdot(\phi f\mathbf{u}) - f\,\mathbf{u}\cdot\nabla\phi$ and the transport
term integrates to zero on a bi-periodic or wall-bounded box. The claim that
makes the header worth having is what follows: the gap between the two is
discretisation error and nothing else, whereas the area alone -- a state
variable -- cannot separate an interface leak from an advection leak. This is
the test of that claim.

A drop of radius $R_0$ is released from a strongly deformed shape,
$r = R_0(1 + \varepsilon\cos 2\theta)$ with $\varepsilon = 0.4$, and damped so
that it relaxes monotonically to its equilibrium circle instead of ringing.
Over the run, $\int\Psi_\sigma\,dt$ is compared against $-\sigma\Delta A$.

Two choices matter, and the obvious alternatives do not work:

* **Large amplitude.** At the $\varepsilon = 0.05$ of
  [oscillation.c](/src/test/oscillation.c) the available area change is
  $2\pi R_0(n^2-1)\varepsilon^2/4$ per unit perimeter, some $10^{-3}$
  relative, which is *below* the drift of `interface_area()` itself: measured
  there, the area passed below the equal-area circle's perimeter, which is
  geometrically impossible and so is pure numerical drift. At
  $\varepsilon = 0.4$ the signal is ~11% and two orders clear of it.
* **Damping.** Undamped, the drop oscillates and $\Delta A$ depends entirely
  on the phase at which the run stops; the comparison then measures the
  stopping time, not the identity. Relaxation to equilibrium makes the
  endpoint unambiguous.

The domain is a quarter of the drop with the default symmetry conditions,
which give $\mathbf{u}\cdot\mathbf{n} = 0$ on every boundary -- the condition
the identity's derivation needs. Both $\Psi_\sigma$ and $A$ are then quarter
quantities and the identity is unaffected.

`interfacial_power_energy.asc` gets one row per resolution,
`N err_area resid`:

* `err_area` compares the final area with the quarter perimeter of the
  equal-area circle, $2\pi R_0\sqrt{1 + \varepsilon^2/2}/4$, known in closed
  form. It is the guard that the drop actually reached equilibrium; without it
  the residual below would be comparing against an arbitrary endpoint.
* `resid` is $|\int\Psi_\sigma dt + \sigma\Delta A|/|\sigma\Delta A|$.

`test_interfacial_power.py` applies the thresholds. Measured, at
$D/\Delta = 12.8$, $25.6$ and $51.2$: `resid` = 0.128, 0.0331, 0.0130, i.e.
observed orders 1.95 then 1.35. The convergence is real but not clean, and
three resolutions cannot separate a genuine order near 1.4 from a floor being
approached near 1%; the checker therefore asserts first order, which is what
the data supports, rather than the second order asserted elsewhere in this
directory.

So the identity holds to about 1% at $D/\Delta = 51$ and improves with
resolution. That number is the useful output: it is the size below which a
non-zero $\Psi_\sigma + \sigma\,dA/dt$ in a production run should be read as
discretisation error rather than as an interface leak.

Only `multigrid` is registered. The three-dimensional analogue -- a spheroid
relaxing viscously to a sphere -- costs far more than the other tests in this
group at a resolution where the residual is meaningful, and 3D coverage of the
routines themselves is provided by the static, routes and AMR tests. That is a
cost decision, not a statement that the identity is two-dimensional. */

#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "tension.h"
#include "acastillo/output_fields/interfacial_power.h"

#define R0 0.1
#define EPS 0.4
#define SIGMA 1.0

int LEVEL;
double work, area0, areaT;
FILE * fp = NULL;

int main() {

  /**
  Damped enough to relax within the time run, with a density ratio so that
  `alpha` varies across the interface. */

  rho1 = 1., rho2 = 1e-2;
  mu1 = 0.05, mu2 = 5e-4;
  f.sigma = SIGMA;
  L0 = 0.5;
  TOLERANCE = 1e-4;

  if (pid() == 0) {
    fp = fopen ("interfacial_power_energy.asc", "w");
    fprintf (fp, "# N err_area resid\n");
  }

  for (LEVEL = 5; LEVEL <= 7; LEVEL++) {
    N = 1 << LEVEL;
    init_grid (N);
    run();
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_interfacial_power_energy.ref`. */
    system ("python3 ../test_interfacial_power.py "
            "--energy interfacial_power_energy.asc 1>&2");
  }
}

event init (i = 0) {
  fraction (f, R0*(1. + EPS*cos (2.*atan2 (y, x))) - sqrt (sq(x) + sq(y)));
  work = 0.;
  area0 = -1.;
}

event projection (i++) {
  /**
  Route 2 is used, but from the same `projection` hook as
  [test_interfacial_power_routes.c](test_interfacial_power_routes.c), before
  `centered.h`'s `correction()`. Route 2 does not *require* that placement --
  it needs neither `a` nor gravity -- but $\Psi_\sigma$ is a functional of
  $\mathbf{u}$, and $\mathbf{u}$ changes within the timestep, so an
  accumulated integral is only meaningful once the sampling point is fixed.
  Sampling after the correction instead changes the result by a factor of
  about two. */

  scalar d[];
  double psi = interfacial_power_curvature (f, d);
  double A = interface_area (f);
  if (area0 < 0.)
    area0 = A;
  areaT = A;
  work += psi*dt;
}

event stop (t = 2.) {

  /**
  Volume is conserved, so the equilibrium circle has
  $R_{eq} = R_0\sqrt{1 + \varepsilon^2/2}$; the domain holds a quarter of its
  perimeter. */

  double area_eq = 2.*pi*R0*sqrt (1. + sq(EPS)/2.)/4.;
  double sdA = -SIGMA*(areaT - area0);

  if (pid() == 0) {
    fprintf (fp, "%d %.17g %.17g\n", N,
             fabs (areaT - area_eq)/area_eq,
             fabs (work - sdA)/fabs (sdA));
    fflush (fp);
  }
  return 1;
}
