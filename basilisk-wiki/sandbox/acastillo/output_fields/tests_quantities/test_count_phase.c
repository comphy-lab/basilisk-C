/**
# Testing `count_phase()` on three spheres

Three well-separated spheres of distinct radii $R_k$ at centres
$\mathbf{c}_k$, in a solid-body rotation about $z$ superposed on a uniform
translation,

$$\mathbf{u} = \mathbf{U}_0 + \Omega\,\hat{z}\times\mathbf{r}$$

give every column of [count_phase.h](../count_phase.h) a closed form:
$V = \tfrac{4}{3}\pi R_k^3$, $\mathbf{b} = \mathbf{c}_k$,
$A = 4\pi R_k^2$, and -- since a sphere is symmetric about its centre, so the
volume-weighted average of a linear field is that field at the centre --
$\langle\mathbf{u}\rangle = \mathbf{U}_0 + \Omega\,\hat{z}\times\mathbf{c}_k$.
Distinct radii let the checker match regions to spheres by volume, `tag()`
numbering them in traversal order.

In 2D the same centres and radii describe three disjoint discs, with
$V = \pi R_k^2$ and $A = 2\pi R_k$ (a perimeter); the checker picks the right
closed form from the row width.

This is the `above = 1` (droplet) branch;
[test_count_phase_bubbles.c](test_count_phase_bubbles.c) covers `above = 0` on
the same geometry.

`count_phase()` writes `count_droplets.asc` itself, so nothing is computed
here: `test_count_phase.py` holds the closed form and the tolerances, and
mirrors the geometry below. */

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
  // Refine near the interfaces, where the VOF reconstruction sets the accuracy
  // of the volume and especially the area.
  double band = 4. * L0 / (1 << LEVEL);
  refine (fabs (phi_spheres (x, y, z)) < band && level < MAXLEVEL);
#endif

  fraction (f, phi_spheres (x, y, z));

  foreach() {
    u.x[] = U0_x - Omega * y;
    u.y[] = U0_y + Omega * x;
#if dimension == 3
    u.z[] = U0_z;
#endif
  }
  boundary ((scalar *){u});

  scalar m[];
  count_phase (m, 1e-6, 1, "count_droplets.asc", 0, 0.);

  /**
  The report goes to stderr, i.e. to the `log` diffed against
  `test_count_phase.ref`. */
  if (pid() == 0)
    system ("python3 ../test_count_phase.py "
            "--droplets count_droplets.asc 1>&2");
}
