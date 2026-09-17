/**
# Testing `interfacial_power_curvature()` on a static drop

A drop of radius $R$ with $f = 1$ inside has a constant potential
$\phi = \sigma\kappa$, so [interfacial_power.h](../interfacial_power.h)'s

$$
\Psi_\sigma = \int_\Omega \mathbf{u}\cdot\phi\nabla f\,dV
  = \sigma\kappa\int_\Omega \mathbf{u}\cdot\nabla f\,dV
  = -\sigma\kappa\oint (\mathbf{u}\cdot\mathbf{n})\,dA
$$

has a closed form for any $\mathbf{u}$, the sign following from $\nabla f$
pointing inwards. Two are used, with the drop clear of the boundaries so that
$\nabla f$ vanishes there:

* $\mathbf{u} = \mathbf{U}$ uniform gives $\Psi_\sigma = 0$, since
  $\oint\mathbf{n}\,dA = 0$ over a closed interface;
* $\mathbf{u} = \alpha\mathbf{x}$ gives $\mathbf{u}\cdot\mathbf{n} = \alpha R$
  and $\Psi_\sigma = -\sigma\kappa\,\alpha R\,A$, which is $-2\pi\sigma\alpha R$
  in 2D ($\kappa = 1/R$, $A = 2\pi R$) and $-8\pi\sigma\alpha R^2$ in 3D
  ($\kappa = 2/R$, $A = 4\pi R^2$).

The drop centre is offset from the origin. On a centred drop the null is
satisfied to roundoff by mirror symmetry alone, whatever the stencil does, and
tests nothing. The offset leaves the second closed form unchanged, since
$\alpha\mathbf{x} = \alpha(\mathbf{x} - \mathbf{x}_c) + \alpha\mathbf{x}_c$ and
the constant part contributes nothing.

The fraction is initialised with [Vofi](#vofi) rather than `fraction()`.
`fraction()` linearises the level set within each cell, leaving an
$O(\Delta/R)$ error in $f$; the height function sums $f$ over a column and the
curvature is a second difference of that divided by $\Delta^2$, so the relative
curvature error settles at $O(1/12)$ *independently of resolution*. Measured on
this drop, $\kappa$ from `fraction()` has an rms error of 5.4% at every
resolution from $D/\Delta = 19$ to $307$, and $\Psi_\sigma$ then plateaus near
$5\times 10^{-3}$ instead of converging. Sub-refining buys a constant factor,
not an order.

`interfacial_power_static.asc` gets one row per resolution,
`N err_uniform err_radial err_field err_linear`, each normalised by the exact
$|\Psi_\sigma|$ of the radial case:

* `err_uniform` and `err_radial` are the two closed forms above, both
  second-order convergent;
* `err_field` checks that the output field `d` integrates to the return value,
  as the header claims;
* `err_linear` checks $\Psi_\sigma(\mathbf{u}_1 + \mathbf{u}_2) =
  \Psi_\sigma(\mathbf{u}_1) + \Psi_\sigma(\mathbf{u}_2)$, which catches a
  `curvature()` call accumulating into $\phi$ or a stale field surviving
  between calls.

Both are roundoff-level. `test_interfacial_power.py` applies the thresholds.

Only the second route is exercised here: `interfacial_power_acceleration()`
reads back `a`, which no solver has filled, and is covered by
[test_interfacial_power_routes.c](test_interfacial_power_routes.c).

## Vofi

[Vofi](https://github.com/VOFTracking/Vofi) computes the volume fraction by
exact quadrature. Its reference phase is where the implicit function is
*negative*, the opposite of `fraction()`, so `drop()` is negative inside and
$f = 1$ there. The current library exports the prefixed names used below; code
written against the pre-2015 API calls them `creal`, `Get_fh` and `Get_cc`. */

#include <vofi.h>
#pragma autolink -L$HOME/local/lib -lvofi

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "acastillo/output_fields/interfacial_power.h"

#define RDROP 0.3
#define SIGMA 1.3
#define ALPHA 0.7

// Uniform velocity of the null case, with no vanishing component.
#define U0_x 0.9
#define U0_y (-0.5)
#define U0_z 0.4

// Offset of the drop centre, breaking the mirror symmetry of the grid.
static double xc = 0.07, yc = -0.043, zc = 0.031;

static double drop (vofi_creal p[dimension])
{
#if dimension == 3
  return sq (p[0] - xc) + sq (p[1] - yc) + sq (p[2] - zc) - sq (RDROP);
#else
  return sq (p[0] - xc) + sq (p[1] - yc) - sq (RDROP);
#endif
}

static void vofi_fraction (scalar c, int Nres)
{
  double fh = vofi_Get_fh (drop, NULL, 1./Nres, dimension, 0);
  foreach() {
    vofi_creal p[3] = {x - Delta/2., y - Delta/2., z - Delta/2.};
    c[] = vofi_Get_cc (drop, p, Delta, fh, dimension);
  }
}

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  f.sigma = SIGMA;

  /**
  $\kappa$ is the sum of the principal curvatures, so $1/R$ for a disc and
  $2/R$ for a sphere; $A$ is a perimeter in 2D and an area in 3D. */

#if dimension == 3
  double kappa = 2./RDROP, area = 4.*pi*sq (RDROP);
  int Nmin = 32, Nmax = 128;
#else
  double kappa = 1./RDROP, area = 2.*pi*RDROP;
  int Nmin = 64, Nmax = 512;
#endif
  double psi_exact = -SIGMA*kappa*ALPHA*RDROP*area;

  FILE * fp = NULL;
  if (pid() == 0) {
    fp = fopen ("interfacial_power_static.asc", "w");
    fprintf (fp, "# N err_uniform err_radial err_field err_linear\n");
  }

  for (int Nres = Nmin; Nres <= Nmax; Nres *= 2) {
    init_grid (Nres);
    vofi_fraction (f, Nres);

    scalar d[];

    foreach() {
      u.x[] = U0_x;
      u.y[] = U0_y;
#if dimension == 3
      u.z[] = U0_z;
#endif
    }
    double psi_uniform = interfacial_power_curvature (f, d);

    foreach() {
      u.x[] = ALPHA*x;
      u.y[] = ALPHA*y;
#if dimension == 3
      u.z[] = ALPHA*z;
#endif
    }
    double psi_radial = interfacial_power_curvature (f, d);

    // The same integral, taken from the field the call leaves behind.
    double psi_field = 0.;
    foreach (reduction(+:psi_field))
      psi_field += dv()*d[];

    // The sum of the two velocity fields, for linearity.
    foreach() {
      u.x[] = U0_x + ALPHA*x;
      u.y[] = U0_y + ALPHA*y;
#if dimension == 3
      u.z[] = U0_z + ALPHA*z;
#endif
    }
    double psi_sum = interfacial_power_curvature (f, d);

    if (pid() == 0) {
      fprintf (fp, "%d %.17g %.17g %.17g %.17g\n", Nres,
               fabs (psi_uniform)/fabs (psi_exact),
               fabs (psi_radial - psi_exact)/fabs (psi_exact),
               fabs (psi_field - psi_radial)/fabs (psi_exact),
               fabs (psi_sum - psi_uniform - psi_radial)/fabs (psi_exact));
      fflush (fp);
    }
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_interfacial_power_static.ref`. */
    system ("python3 ../test_interfacial_power.py "
            "--static interfacial_power_static.asc 1>&2");
  }
}
