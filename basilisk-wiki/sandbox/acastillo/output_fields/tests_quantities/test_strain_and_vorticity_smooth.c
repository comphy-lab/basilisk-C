/**
# Testing `strain_and_vorticity_sq()` on a smooth field

Where
[test_strain_and_vorticity_affine.c](test_strain_and_vorticity_affine.c)
removes the truncation error to isolate the algebra, this checks the
truncation error itself: on a smooth field the centred differences in
[strain_and_vorticity.h](../strain_and_vorticity.h) are second-order, so the
error must fall as $\Delta^2$.

Only $\Omega^2$ is measured here. Its $S^2$ is bit-for-bit
[strain_rate_sq()](../strain_rate.h)'s -- the affine test asserts exactly
that -- so the convergence of $S^2$ is already covered by
[test_strain_rate_smooth.c](test_strain_rate_smooth.c), and repeating it would
freeze a second reference for the same numbers.

The fields are those of `test_strain_rate_smooth.c`. In 2D, from
$\psi = \cos(k_1 x)\cos(k_2 y)$ with
$\mathbf{u} = (\partial_y\psi, -\partial_x\psi)$, the only vorticity component
is $\omega_z = -\nabla^2\psi = (k_1^2 + k_2^2)\cos(k_1x)\cos(k_2y)$. In 3D the
ABC flow is a Beltrami field, $\boldsymbol\omega = k\mathbf{u}$. Either way

$$\Omega_{ij}\Omega_{ij} = \tfrac{1}{2}|\boldsymbol\omega|^2$$

which is what the exact value below evaluates -- an independent route to
$\Omega^2$, not a transcription of the code's expression.

The 3D table comes out equal to
[test_strain_rate_smooth.c](test_strain_rate_smooth.c)'s, which is expected
rather than a stray copy: the ABC flow has
$S_{ij}S_{ij} + \Omega_{ij}\Omega_{ij} = k^2(A^2 + B^2 + C^2)$, a constant, so
the two errors differ only in sign and share a normalisation. The fields
themselves are far apart -- $S^2$ and $\Omega^2$ differ by up to 529 here --
so the test still fails loudly if $\Omega^2$ is wrong. The 2D table does not
coincide.

Cells are inset one layer from the boundary, since the stencil would otherwise
reach a ghost cell -- which is also why no periodicity is imposed.
`strain_and_vorticity_smooth.asc` gets one row per resolution,
`N err_max err_l2`, both normalised by $\max|\Omega^2|$;
`test_strain_and_vorticity.py` turns them into convergence orders. */

#include "run.h"
#include "acastillo/output_fields/strain_and_vorticity.h"

double k1, k2, kabc;
double Aabc = 1.0, Babc = 2.0, Cabc = 3.0;

double O2_smooth (double x, double y, double z) {
#if dimension == 3
  // Beltrami: omega = kabc*u, so O2 = |omega|^2/2.
  double ux = Aabc*sin(kabc*z) + Cabc*cos(kabc*y);
  double uy = Babc*sin(kabc*x) + Aabc*cos(kabc*z);
  double uz = Cabc*sin(kabc*y) + Babc*cos(kabc*x);
  return 0.5*sq(kabc)*(sq(ux) + sq(uy) + sq(uz));
#else
  double omega = (sq(k1) + sq(k2))*cos(k1*x)*cos(k2*y);
  return 0.5*sq(omega);
#endif
}

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;

  k1   = 2.*pi/L0;
  k2   = 4.*pi/L0;
  kabc = 2.*pi/L0;

  FILE * fp = NULL;
  if (pid() == 0) {
    fp = fopen ("strain_and_vorticity_smooth.asc", "w");
    fprintf (fp, "# N err_max err_l2\n");
  }

  // Three resolutions give two independent convergence rates. The coarsest 3D
  // grid still resolves the single wavelength with 8 cells.
#if dimension == 3
  int Nmin = 8, Nmax = 32;
#else
  int Nmin = 32, Nmax = 128;
#endif

  for (int Nres = Nmin; Nres <= Nmax; Nres *= 2) {
    init_grid (Nres);

    vector u[];

    foreach() {
#if dimension == 3
      u.x[] = Aabc*sin(kabc*z) + Cabc*cos(kabc*y);
      u.y[] = Babc*sin(kabc*x) + Aabc*cos(kabc*z);
      u.z[] = Cabc*sin(kabc*y) + Babc*cos(kabc*x);
#else
      u.x[] = -k2*cos(k1*x)*sin(k2*y);
      u.y[] =  k1*sin(k1*x)*cos(k2*y);
#endif
    }

    double err_max = 0., err_sq = 0., o2_max = 0., vol = 0.;
    foreach (reduction(max:err_max) reduction(+:err_sq)
             reduction(max:o2_max) reduction(+:vol)) {
      bool interior = fabs(x) < L0/2. - 1.1*Delta && fabs(y) < L0/2. - 1.1*Delta;
#if dimension == 3
      interior = interior && fabs(z) < L0/2. - 1.1*Delta;
#endif
      if (interior) {
        double S2, O2;
        strain_and_vorticity_sq (point, u, &S2, &O2);
        double o2_exact = O2_smooth (x, y, z);
        double e = fabs (O2 - o2_exact);
        if (e > err_max)
          err_max = e;
        if (fabs(o2_exact) > o2_max)
          o2_max = fabs(o2_exact);
        err_sq += dv()*sq(e);
        vol += dv();
      }
    }

    if (pid() == 0) {
      fprintf (fp, "%d %.17g %.17g\n", Nres, err_max/o2_max,
               sqrt(err_sq/vol)/o2_max);
      fflush (fp);
    }
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_strain_and_vorticity_smooth.ref`. */
    system ("python3 ../test_strain_and_vorticity.py "
            "--smooth strain_and_vorticity_smooth.asc 1>&2");
  }
}
