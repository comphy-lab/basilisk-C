
void WARN_m_dot(double m_dot_sum_pre, double m_dot_sum_post) {
  double m_dot_diff =
    fabs(m_dot_sum_post - m_dot_sum_pre) / max(fabs(m_dot_sum_pre), F_ERR);
  static int WARN = 1;
  if (m_dot_diff > 1.0e-3) {
    if (WARN == 1) {
      WARN++;
      fprintf(stderr, "Warning: m_dot_diff>1.0e-3 %g\n", m_dot_diff);
      fprintf(stderr, "m_dot_sum_post= %g m_dot_sum_pre=%g\n",
              m_dot_sum_post, m_dot_sum_pre);
    }
  }
}

#define WARNING_FLAG(err)                                                \
    if (FLAG == 0) {                                                     \
        static int ONE_WARN = 1;                                         \
        if (ONE_WARN == 1) {                                             \
            ONE_WARN = 0;                                                \
            fprintf(stderr,                                              \
                "%d Warning: no close interface cell found to m_dot in " \
                "pure cell "                                             \
                "directon\n",                                            \
                err);                                                    \
            fprintf(stderr, "x= %g y=%g z=%g\n", x, y, z);               \
        }                                                                \
    }

void donate_mdot_func(scalar m_dot_hat, vector n_vec, scalar interfacial_point, bool VAP) {
  // translate m_dot location to nearest pure gas cell
  //////////////////////////////////
  // IMPORTANT BC call because these will be accesses using noauto loop with
  // foreach_neighbor
  m_dot.dirty = true;
  m_dot_hat.dirty = true;
  f.dirty = true;
  interfacial_point.dirty = true;
  boundary((scalar *){m_dot, m_dot_hat, f, interfacial_point});
  //////////////////////////////////////////////////

#if _MPI
  MPI_DIST_INIT;
#endif

  double m_dot_sum_pre = 0.0;
  foreach (reduction(+:m_dot_sum_pre)) {
    m_dot_sum_pre += m_dot[];
  }

  foreach (noauto) {
    if (INTERFACIAL && fabs(m_dot[]) > F_ERR) {
      // if pure vapor then you dont need to distrbute m_dot
      // if (NOT_INTERFACIAL && is_VAPOR) {
      //    printf("now not interfacial but m_dot is real");
      //}else{
      coord cell_center = {x, y, z};
      coord n_interface = {0., 0., 0.};
      foreach_dimension() { n_interface.x = n_vec.x[]; }

      double eps_val[cells_G];
      initialize_stencil(eps_val);
      //
      int FLAG = 0;
      int index = -1;
      // this is why we need noauto
      foreach_neighbor() {
        index++;
        coord cell_current = {x, y, z};
        // not interface cell - pure, vapor and not center of stencil
        if NOT_BOUNDARY {
          //if (NOT_INTERFACIAL && is_VAPOR) {
          if (NOT_INTERFACIAL) {
            if (VAP ? is_VAPOR : is_LIQUID){
              double dist = 0.0;

              //get_collin_factor(true, true, n_interface,
              double eps_col =
                  get_collin_factor(VAP, true, n_interface,
                                    cell_center, cell_current, &dist);
              eps_val[index] = eps_col;
              FLAG = 1;
            }
          }
        }
      }

      WARNING_FLAG(2);
      //
      // normalize the distance values for distribution
      normalize_stencil(eps_val);
      //
      double m_dot_current = m_dot[];

      index = -1;
      foreach_neighbor() {
          index++;
          if (eps_val[index] > 1.0e-6) {
            double val_mdot = eps_val[index] * m_dot_current;
            //
#if _MPI
            MPI_DIST_SUM(m_dot_hat, val_mdot, 53);
#else
            m_dot_hat[] += val_mdot;
#endif  // _MPI
          }
      }
      // set interface m_dot to zero
      // commented because now using m_dot_hat
      //m_dot[] = 0.0;
      //}
    }
  }
  
#if _MPI
  MPI_DIST_FINISH_SUM(m_dot_hat);
#endif  // _MPI
  //
  m_dot_hat.dirty = true;
  boundary((scalar *){m_dot_hat});

  double m_dot_sum_post = 0.0;
  foreach (reduction(+:m_dot_sum_post)) {
    m_dot_sum_post += m_dot_hat[];
  }

#if AXI
  foreach () {
    if (fabs(m_dot[]) > F_ERR) {
      //m_dot[] *= y;
      m_dot_hat[] *= y;
    }
  }
#endif

  WARN_m_dot(m_dot_sum_pre, m_dot_sum_post);

  return;
}

void compute_mdot(vector n_vec, scalar interfacial_point) {
  foreach () {
    m_dot[] = 0.0;
    j_dot[] = 0.0;
  }
  foreach () {
    if INTERFACIAL {
      coord n_interface = {0., 0., 0.};
      coord p = {0., 0., 0.};
      foreach_dimension() { n_interface.x = n_vec.x[]; }
      double alpha = plane_alpha(f[], n_interface);
      double a_interface = plane_area_center(n_interface, alpha, &p);
#if AXI
      if (y > 1.0e-12) {
        a_interface = a_interface * (y + Delta * p.y) / y;
      }
#endif  // AXI
  // interfacial area density [m^2/m^3]
      double iad = a_interface / Delta;
      //
      double gradient_sum =
          lambda_L * dT_interface_L[] - lambda_V * dT_interface_V[];
      double j_dot_val = gradient_sum / hlg;
      j_dot[] = j_dot_val;
      m_dot[] = j_dot_val * iad;
      //}
    } else {
      m_dot[] = 0.0;
      j_dot[] = 0.0;
    }
  }
}

void distribute_mdot(vector n_vec, scalar interfacial_point) {

  scalar m_dot_hat[];
  foreach(){
    m_dot_hat[] = 0.0;
  }

  bool VAP = true;
  donate_mdot_func(m_dot_hat, n_vec, interfacial_point, VAP);
  VAP = false;
  donate_mdot_func(m_dot_hat, n_vec, interfacial_point, VAP);

  foreach(){
    m_dot[] = m_dot_hat[];
  }
}
