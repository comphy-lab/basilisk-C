double global_fm_min = HUGE;
double timestep(const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  double fm_min = HUGE;
  foreach_face(reduction(min : dtmax)) if (u.x[] != 0.)
  {
    double dt = Delta / fabs(u.x[]);
    fm_min = fmin(fm_min, fm.x[]);
#if EMBED
    assert(fm.x[]);
    dt *= fm.x[];
#else
    dt *= cm[];
#endif
    if (dt < dtmax)
      dtmax = dt;
  }
  dtmax *= CFL;

  global_fm_min = fmin(fm_min, global_fm_min);
  fprintf(stderr, "fm_min=%.6f, glob_fm_min=%.6f\n", fm_min, global_fm_min);
  fflush(stderr);

  /**
  We rest *previous* between successive runs. */

  if (t == 0.)
    previous = 0.;

  if (dtmax > previous)
    dtmax = (previous + 0.1 * dtmax) / 1.1;
  previous = dtmax;
  return dtmax;
}