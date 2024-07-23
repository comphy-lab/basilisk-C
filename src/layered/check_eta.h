/**
# Consistency check on free-surface evolution

This is optional and only applies to the [time-implicit](implicit.h)
layered solver.

This checks that the free-surface position $\eta$ obtained from the
solution of the implicit time-integration matches the free-surface
position $\eta_p$ obtained by integration of the corresponding
(barotropic) fluxes i.e.
$$
\partial_t\eta_p + \sum_l \nabla\cdot(\mathbf{hu})_l = 0
$$
This should be true to within the tolerance of the implicit solver.

The difference is stored in the `deta` field. */

scalar deta[];

event update_eta (i++)
{
  foreach() {
    double eta_p = zb[];
    foreach_layer()
      eta_p += h[];
    deta[] = eta_p - eta[];
    if (fabs(deta[]) > 1.1*TOLERANCE)
      fprintf (stderr, "warning: fabs(etap - eta[]) = %g > 1.1*TOLERANCE = %g "
	       "at %g,%g,%g\n",
	       fabs(deta[]), 1.1*TOLERANCE, x, y, t);
  }
}
