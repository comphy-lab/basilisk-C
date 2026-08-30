/**
# Testing `strain_rate_sq()` on a smooth field

Where [test_strain_rate_affine.c](test_strain_rate_affine.c) removes the
truncation error to isolate the algebra, this checks the truncation error
itself: on a smooth field the centred differences in
[strain_rate.h](../strain_rate.h) are second-order, so the error must fall as
$\Delta^2$.

2D, from $\psi = \cos(k_1 x)\cos(k_2 y)$ with
$\mathbf{u} = (\partial_y\psi, -\partial_x\psi)$:

$$S_{xx} = -S_{yy} = k_1k_2\sin(k_1x)\sin(k_2y), \qquad
  S_{xy} = \tfrac{1}{2}(k_1^2-k_2^2)\cos(k_1x)\cos(k_2y)$$

$k_1 \neq k_2$ is required: otherwise $S_{xy} \equiv 0$ and the off-diagonal
strain goes untested. 3D uses the ABC flow, which has
$S_{xx} = S_{yy} = S_{zz} = 0$ and all three off-diagonal terms non-zero;
distinct $A, B, C$ keep them from coinciding.

Cells are inset one layer from the boundary, since the stencil would otherwise
reach a ghost cell -- which is also why no periodicity is imposed.
`strain_rate_smooth.asc` gets one row per resolution, `N err_max err_l2`, both
normalised by $\max|S^2|$; `test_strain_rate.py` turns them into convergence
orders. */

#include "run.h"
#include "acastillo/output_fields/strain_rate.h"

double k1, k2, kabc;
double Aabc = 1.0, Babc = 2.0, Cabc = 3.0;

double S2_smooth (double x, double y, double z) {
#if dimension == 3
  double Sxy = 0.5*kabc*(Babc*cos(kabc*x) - Cabc*sin(kabc*y));
  double Sxz = 0.5*kabc*(Aabc*cos(kabc*z) - Babc*sin(kabc*x));
  double Syz = 0.5*kabc*(Cabc*cos(kabc*y) - Aabc*sin(kabc*z));
  return 2.*(sq(Sxy) + sq(Sxz) + sq(Syz));
#else
  double Sxx = k1*k2*sin(k1*x)*sin(k2*y);
  double Sxy = 0.5*(sq(k1) - sq(k2))*cos(k1*x)*cos(k2*y);
  return 2.*(sq(Sxx) + sq(Sxy));
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
    fp = fopen ("strain_rate_smooth.asc", "w");
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

    double err_max = 0., err_sq = 0., s2_max = 0., vol = 0.;
    foreach (reduction(max:err_max) reduction(+:err_sq)
             reduction(max:s2_max) reduction(+:vol)) {
      bool interior = fabs(x) < L0/2. - 1.1*Delta && fabs(y) < L0/2. - 1.1*Delta;
#if dimension == 3
      interior = interior && fabs(z) < L0/2. - 1.1*Delta;
#endif
      if (interior) {
        double s2_exact = S2_smooth (x, y, z);
        double e = fabs (strain_rate_sq (point, u) - s2_exact);
        if (e > err_max)
          err_max = e;
        if (fabs(s2_exact) > s2_max)
          s2_max = fabs(s2_exact);
        err_sq += dv()*sq(e);
        vol += dv();
      }
    }

    if (pid() == 0) {
      fprintf (fp, "%d %.17g %.17g\n", Nres, err_max/s2_max,
               sqrt(err_sq/vol)/s2_max);
      fflush (fp);
    }
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_strain_rate_smooth.ref`. */
    system ("python3 ../test_strain_rate.py "
            "--smooth strain_rate_smooth.asc 1>&2");
  }
}
