/**
# Testing `strain_and_vorticity_sq()` on an affine field

On $u_i = A_{ij}x_j$ with $\mathrm{tr}(A) = 0$, the gradient is constant, so
the centred differences in
[strain_and_vorticity.h](../strain_and_vorticity.h) are exact and both

$$S^2 = S_{ij}S_{ij}, \qquad S_{ij} = \tfrac{1}{2}(A_{ij} + A_{ji})$$
$$\Omega^2 = \Omega_{ij}\Omega_{ij}, \qquad
  \Omega_{ij} = \tfrac{1}{2}(A_{ij} - A_{ji})$$

are single constants. Any error is then algebraic, not truncation. The same
matrices as [test_strain_rate_affine.c](test_strain_rate_affine.c) are used,
so the shared $S^2$ is compared against the same reference value; they are
chosen with $A_{ij} \neq A_{ji}$ throughout, or $\Omega^2$ would vanish and go
untested.

The third column checks the claim in the header that this $S^2$ is the same
quantity as [strain_rate_sq()](../strain_rate.h)'s, the two living side by
side. They are not bit-for-bit: the off-diagonal terms are summed as
$2S_{xy}^2$ in one and $S_{xy}^2 + S_{yx}^2$ in the other, so the two round
differently and the check is against roundoff, not against zero. What it
catches is the pair drifting apart algebraically.

Cells are inset one layer from the boundary, since the stencil would otherwise
reach a ghost cell. `strain_and_vorticity_affine.asc` gets one row per
resolution, `N err_S2 err_O2 err_cons`, the first two the max relative
deviation from the constants; `test_strain_and_vorticity.py` applies the
thresholds. */

#include "run.h"
#include "acastillo/output_fields/strain_and_vorticity.h"
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

  // Contracted from A_aff, so they cannot drift from the coefficients.
  double s2_exact = 0., o2_exact = 0.;
  for (int i = 0; i < dimension; i++)
    for (int j = 0; j < dimension; j++) {
      s2_exact += sq (0.5*(A_aff[i][j] + A_aff[j][i]));
      o2_exact += sq (0.5*(A_aff[i][j] - A_aff[j][i]));
    }

  FILE * fp = NULL;
  if (pid() == 0) {
    fp = fopen ("strain_and_vorticity_affine.asc", "w");
    fprintf (fp, "# N err_S2 err_O2 err_cons\n");
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

    double err_s2 = 0., err_o2 = 0., err_cons = 0.;
    foreach (reduction(max:err_s2) reduction(max:err_o2)
             reduction(max:err_cons)) {
      bool interior = fabs(x) < L0/2. - 1.1*Delta && fabs(y) < L0/2. - 1.1*Delta;
#if dimension == 3
      interior = interior && fabs(z) < L0/2. - 1.1*Delta;
#endif
      if (interior) {
        double S2, O2;
        strain_and_vorticity_sq (point, u, &S2, &O2);
        double e = fabs (S2 - s2_exact);
        if (e > err_s2)
          err_s2 = e;
        e = fabs (O2 - o2_exact);
        if (e > err_o2)
          err_o2 = e;
        e = fabs (S2 - strain_rate_sq (point, u));
        if (e > err_cons)
          err_cons = e;
      }
    }

    if (pid() == 0) {
      fprintf (fp, "%d %.17g %.17g %.17g\n", Nres, err_s2/s2_exact,
               err_o2/o2_exact, err_cons/s2_exact);
      fflush (fp);
    }
  }

  if (pid() == 0) {
    fclose (fp);

    /**
    The report goes to stderr, i.e. to the `log` diffed against
    `test_strain_and_vorticity_affine.ref`. */
    system ("python3 ../test_strain_and_vorticity.py "
            "--affine strain_and_vorticity_affine.asc 1>&2");
  }
}
