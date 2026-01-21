#pragma autolink -lgsl -lgslcblas
#include <gsl/gsl_poly.h>

// Adaptive timestep with quadratic CFL: a0*dt^2 + |u|*dt - Delta = 0
double timestep_force (const face vector u, double dtmax, double a0)
{
  static double previous2 = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt1, dt2;
      // Solve quadratic for timestep accounting for acceleration
      gsl_poly_solve_quadratic(a0, fabs(u.x[]), Delta, &dt1, &dt2);
      double dt = max(dt1, dt2);

      dt *= cm[];  // Metric correction

      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  // Smooth timestep growth (max 10% increase per step)
  if (dtmax > previous2)
    dtmax = (previous2 + 0.1*dtmax)/1.1;
  previous2 = dtmax;
  return dtmax;
}