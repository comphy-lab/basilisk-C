/**
# Testing `count_phase()` on three bubbles

Same three spheres and the same velocity field as
[test_count_phase.c](test_count_phase.c), but with `f` inverted so the spheres
are the low-`f` phase, counted through the `above = 0` branch. Every closed
form is unchanged, so the two tests must agree to their tolerances.

That agreement is the point. In the `above = 0` branch `count_phase()` forms
`fval = 1 - f[]` and passes it to `plane_alpha()` together with a normal taken
from `f` rather than from `fval`. The plane it reconstructs is therefore not
the interface of the counted phase but its mirror image through the cell
centre -- which has the same area, so `area` should still come out right. This
test is what makes that a verified claim rather than an assumption.

`test_count_phase.py` holds the closed form and mirrors the geometry. 3D only
-- `count_phase()` uses `u.z`. */

#define LEVEL 5
#define MAXLEVEL 7

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "acastillo/output_fields/count_phase.h"

double R_sph[3] = {0.18, 0.14, 0.10};
coord c_sph[3] = {{-0.25, -0.25, -0.25},
                  { 0.25,  0.25, -0.25},
                  { 0.20, -0.28,  0.28}};

double U0_x = 0.1, U0_y = 0.2, U0_z = 0.3;
double Omega = 1.5;

// Level set of the union: positive inside any sphere.
double phi_spheres (double x, double y, double z) {
  double phi = -HUGE;
  for (int k = 0; k < 3; k++) {
    double d2 = sq(x - c_sph[k].x) + sq(y - c_sph[k].y);
#if dimension == 3
    d2 += sq(z - c_sph[k].z);
#endif
    double d = sqrt (d2);
    if (R_sph[k] - d > phi)
      phi = R_sph[k] - d;
  }
  return phi;
}

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid (N);

#if TREE
  double band = 4. * L0 / (1 << LEVEL);
  refine (fabs (phi_spheres (x, y, z)) < band && level < MAXLEVEL);
#endif

  // Inverted: f = 0 inside the spheres, 1 in the ambient fluid.
  fraction (f, -phi_spheres (x, y, z));

  foreach() {
    u.x[] = U0_x - Omega * y;
    u.y[] = U0_y + Omega * x;
#if dimension == 3
    u.z[] = U0_z;
#endif
  }
  boundary ((scalar *){u});

  scalar m[];
  count_phase (m, 1e-6, 0, "count_bubbles.asc", 0, 0.);

  /**
  The report goes to stderr, i.e. to the `log` diffed against
  `test_count_phase_bubbles.ref`. */
  if (pid() == 0)
    system ("python3 ../test_count_phase.py "
            "--bubbles count_bubbles.asc 1>&2");
}
