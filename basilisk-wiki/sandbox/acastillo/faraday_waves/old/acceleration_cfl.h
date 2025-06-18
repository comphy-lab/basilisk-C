#pragma autolink -lgsl -lgslcblas
#include <gsl/gsl_poly.h>
double timestep_force (const face vector u, double dtmax, double a0)
{
  static double previous2 = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt1, dt2;
      gsl_poly_solve_quadratic(a0, fabs(u.x[]), Delta, &dt1, &dt2);
      double dt = max(dt1, dt2);

      dt *= cm[];

      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous2)
    dtmax = (previous2 + 0.1*dtmax)/1.1;
  previous2 = dtmax;
  return dtmax;
}

event set_dtmax (i++,last) {
  dtmax = timestep_force (u, DT, force.Gn*force.ramp);
}
