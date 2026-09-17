/**
# Testing the two routes of `interfacial_power.h` on a tree

[test_interfacial_power_routes.c](test_interfacial_power_routes.c) compares the
two routes on a uniform grid, where they agree to roundoff. On a tree they are
*not* equivalent: `iforce.h` swaps `f`'s prolongation to the pressure's before
differencing, to keep the interfacial force well-balanced against the pressure
gradient at refinement boundaries, and restores it afterwards; route 2 rebuilds
$\phi\nabla f$ without that swap. The two therefore differ wherever the
interface sits on a level jump.

[interfacial_power.h](../interfacial_power.h) claims the difference vanishes
while the interface is held at a single resolution. That is the assumption a
caller makes when it refines a band around the interface and reads
$\Psi_\sigma$ from the cheap route, so it is worth pinning from both sides. The
same drop is run twice on the same tree grid:

* **contained** -- a band around the initial interface is refined to
  `maxlevel`, wide enough that the drop stays inside it for the whole run. The
  interface never meets a level jump and the two routes must agree to roundoff,
  as on a uniform grid.
* **straddling** -- the half-domain $x < 0$ is refined instead, putting a level
  jump straight through the drop. The two routes must now differ by a margin
  far above roundoff.

Without the second case the first is vacuous: a comparison that returns zero
because nothing is being differenced would pass it. The gap between the two
columns is the result.

`interfacial_power_amr.asc` gets one row per timestep,
`mode i err_routes psi`, with `mode` 0 for contained and 1 for straddling;
`test_interfacial_power.py` applies the thresholds. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "acastillo/output_fields/interfacial_power.h"

#define RDROP 0.25
#define SIGMA 1.0

#if dimension == 3
# define MINLEVEL 4
# define MAXLEVEL 6
#else
# define MINLEVEL 5
# define MAXLEVEL 7
#endif

// 0: refined band following the interface. 1: level jump through the drop.
int mode;

FILE * fp = NULL;

// Oblate: at rest this is not an equilibrium shape, so the drop oscillates.
double phi_drop (double x, double y, double z) {
  return sq (RDROP) - sq (x/1.3) - sq (1.3*y) - sq (z);
}

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;

  rho1 = 1.;
  rho2 = 0.1;
  mu1 = mu2 = 0.;
  f.sigma = SIGMA;

  if (pid() == 0) {
    fp = fopen ("interfacial_power_amr.asc", "w");
    fprintf (fp, "# mode i err_routes psi\n");
  }

  for (mode = 0; mode <= 1; mode++) {
    init_grid (1 << MINLEVEL);
    run();
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_interfacial_power_amr.ref`. */
    system ("python3 ../test_interfacial_power.py "
            "--amr interfacial_power_amr.asc 1>&2");
  }
}

event init (i = 0) {

  /**
  The grid is built before the fraction, so that `fraction()` resolves the
  interface at the level it will be carried at. The band is eight cells of
  `maxlevel` on either side of the interface, which the drop does not leave in
  the timesteps run here. */

  if (mode == 0) {
    double band = 8.*L0/(1 << MAXLEVEL);
    refine (fabs (phi_drop (x, y, z)) < band*2.*RDROP && level < MAXLEVEL);
  }
  else
    refine (x < 0. && level < MAXLEVEL);

  fraction (f, phi_drop (x, y, z));
}

event projection (i++) {
  coord ag = {0.};
  scalar d[];
  double psi_a = interfacial_power_acceleration (ag);
  double psi_c = interfacial_power_curvature (f, d);

  double scale = fabs (psi_a) > 0. ? fabs (psi_a) : 1.;
  if (pid() == 0) {
    fprintf (fp, "%d %d %.17g %.17g\n", mode, i,
             fabs (psi_a - psi_c)/scale, psi_a);
    fflush (fp);
  }
}

event stop (i = 20) {
  return 1;
}
