
double timestep_given_vel(double vel_abs, double dtmax, double CFL_) {
  vel_abs = max(vel_abs, F_ERR);
  double delta_min = get_delta_min();
  double dt_local = CFL_ * delta_min / vel_abs;
  dtmax = min(dt_local, dtmax);
  return dtmax;
}

event set_dtmax(i++, last) dtmax = DT;

event stability(i++, last) {
  double vel_shift_max = 0.;
  foreach (reduction(max:vel_shift_max)) {
    if (fabs(j_dot[]) > F_ERR) {
      // j_dot ~ mdot*Delta
      double vel_shift = fabs(j_dot[] / rho_L);
      vel_shift_max = max(vel_shift_max, vel_shift);
    }
  }

  double j_dot_max = 0.;
  foreach (reduction(max:j_dot_max)) {
    // j_dot ~ mdot*Delta
    j_dot_max = max(j_dot_max, fabs(j_dot[]));
  }

  double CFL_shift = 0.1;
  dtmax = timestep_given_vel(vel_shift_max, dtmax, CFL_shift);

  double vel_stef_max = 0.;
  foreach (reduction(max:vel_stef_max)) {
    if (fabs(j_dot[]) > F_ERR) {
      // j_dot ~ mdot*Delta
      // double vel_stef = fabs(evap[]*Delta);
      double vel_stef = fabs(evap_const * j_dot[]);
      vel_stef_max = max(vel_stef_max, vel_stef);
    }
  }
  dtmax = timestep_given_vel(vel_stef_max, dtmax, CFL);
}

