/**
# Testing `strain_rate_sq()` on an affine field

On $u_i = A_{ij}x_j$ with $\mathrm{tr}(A) = 0$, the gradient is constant, so
the centred differences in [strain_rate.h](../strain_rate.h) are exact and

$$S^2 = S_{ij}S_{ij}, \qquad S_{ij} = \tfrac{1}{2}(A_{ij} + A_{ji})$$

is a single constant. Any error is then algebraic, not truncation. The trace is
zero because `strain_rate_sq()` does not remove it. Convergence of the
differences themselves is covered by
[test_strain_rate_smooth.c](test_strain_rate_smooth.c).

Cells are inset one layer from the boundary, since the stencil would otherwise
reach a ghost cell. `strain_rate_affine.asc` gets one row per resolution,
`N err_rel`, the max relative deviation from the constant;
`test_strain_rate.py` applies the threshold. */

#include "run.h"
#include "acastillo/output_fields/strain_rate.h"

#if dimension == 3
double A_aff[3][3] = {{ 0.7,  1.3,  0.5},
                      { 0.4, -0.3,  0.9},
                      { 0.2,  0.6, -0.4}};
#else
double A_aff[2][2] = {{ 0.7,  1.3},
                      { 0.4, -0.7}};
#endif

int main() {

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;

  // Contracted from A_aff, so it cannot drift from the coefficients.
  double s2_exact = 0.;
  for (int i = 0; i < dimension; i++)
    for (int j = 0; j < dimension; j++)
      s2_exact += sq (0.5*(A_aff[i][j] + A_aff[j][i]));

  FILE * fp = NULL;
  if (pid() == 0) {
    fp = fopen ("strain_rate_affine.asc", "w");
    fprintf (fp, "# N err_rel\n");
  }

  // Three resolutions, so a flat error is distinguishable from a small one.
#if dimension == 3
  int Nmin = 8, Nmax = 32;
#else
  int Nmin = 32, Nmax = 128;
#endif

  for (int Nres = Nmin; Nres <= Nmax; Nres *= 2) {
    init_grid (Nres);

    vector u[];

    // Component by component: the row of A_aff differs per component, which
    // foreach_dimension() cannot express.
    foreach() {
#if dimension == 3
      double px = x, py = y, pz = z;
      u.x[] = A_aff[0][0]*px + A_aff[0][1]*py + A_aff[0][2]*pz;
      u.y[] = A_aff[1][0]*px + A_aff[1][1]*py + A_aff[1][2]*pz;
      u.z[] = A_aff[2][0]*px + A_aff[2][1]*py + A_aff[2][2]*pz;
#else
      double px = x, py = y;
      u.x[] = A_aff[0][0]*px + A_aff[0][1]*py;
      u.y[] = A_aff[1][0]*px + A_aff[1][1]*py;
#endif
    }

    double err_max = 0.;
    foreach (reduction(max:err_max)) {
      bool interior = fabs(x) < L0/2. - 1.1*Delta && fabs(y) < L0/2. - 1.1*Delta;
#if dimension == 3
      interior = interior && fabs(z) < L0/2. - 1.1*Delta;
#endif
      if (interior) {
        double e = fabs (strain_rate_sq (point, u) - s2_exact);
        if (e > err_max)
          err_max = e;
      }
    }

    if (pid() == 0) {
      fprintf (fp, "%d %.17g\n", Nres, err_max/s2_exact);
      fflush (fp);
    }
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_strain_rate_affine.ref`. */
    system ("python3 ../test_strain_rate.py "
            "--affine strain_rate_affine.asc 1>&2");
  }
}
