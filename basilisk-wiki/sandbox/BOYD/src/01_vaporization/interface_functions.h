#define vecdist2(a, b) \
    (sq((a).x - (b).x) + sq((a).y - (b).y) + sq((a).z - (b).z))
#define vecdist_1(a) (sqrt(sq((a).x) + sq((a).y) + sq((a).z)))
#define vecdot(a, b) ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define vecdiff(a, b) ((coord){(a).x - (b).x, (a).y - (b).y, (a).z - (b).z})

#if dimension == 1
#define norm_coord(v) fabs(v.x)
#elif dimension == 2
#define norm_coord(v) (sqrt(sq(v.x) + sq(v.y)))
#else  // dimension == 3
#define norm_coord(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#endif  // dimension

void compute_normal(scalar f, vector normal_vector) {
  foreach () {
    coord n = interface_normal(point, f);
    foreach_dimension() { normal_vector.x[] = n.x; }
  }
}

void initialize_stencil(double *stencil_) {
  for (int ii = 0; ii < cells_G; ii++) {
    stencil_[ii] = 0.0;
  }
  return;
}

void normalize_stencil(double *stencil_) {
  double sum = 0.0;
  for (int ii = 0; ii < cells_G; ii++) {
    sum += stencil_[ii];
  }
  if (fabs(sum) < 1.0e-6) {
    sum = 1.0e-6;
  }
  // use sum to normalize
  for (int ii = 0; ii < cells_G; ii++) {
    stencil_[ii] /= sum;
  }
  return;
}

void scale_distance(coord *n) {
  double delta_min = get_delta_min();
  foreach_dimension() { n->x /= delta_min; }
}

double get_dist_plus_collin(bool VAP, bool DIRECTION, coord n_interface,
                            coord cell_center, coord cell_current, double *dist,
                            double *dist_scaled, double *dot_val,
                            double *dot_val_scaled) {
  if (!VAP) {
    // looking in the liq direction
    foreach_dimension() { n_interface.x *= -1.0; }
  }
  //
  // this also determines the sign of dirrection of collin
  // coord dist_vec_scaled = vecdiff(cell_current, cell_center);
  coord dist_vec_scaled = vecdiff(cell_current, cell_center);
  // pre scaling
  *dist = vecdist_1(dist_vec_scaled);
  *dot_val = vecdot(dist_vec_scaled, n_interface);
  // needed to be scaled outherwise bot (sq(dot_val will be very small
  scale_distance(&dist_vec_scaled);
  *dist_scaled = vecdist_1(dist_vec_scaled);
  *dot_val_scaled = vecdot(dist_vec_scaled, n_interface);

  // Malan Eq (24)
  *dot_val_scaled = vecdot(dist_vec_scaled, n_interface);
  // double collin = *dot_val_scaled;
  // new
  coord dist_vec_norm =
    dist_vec_scaled;  // vecdiff(cell_center, cell_current);
  if (*dist > F_ERR) {
    normalize(&dist_vec_norm);
  }
  double collin = vecdot(dist_vec_norm, n_interface);
  if (DIRECTION) {
    collin = max(collin, 0.0);
  }
  return collin;
}

// function returns collin / sq(distance)
double get_collin_factor(bool VAP, bool DIRECTION, coord n_interface,
                         coord cell_center, coord cell_current, double *dist) {
    //
  double dot_val = 0.0;
  double dot_val_scaled = 0.0;
  double dist_scaled = 0.0;
  double collin = get_dist_plus_collin(
      VAP, DIRECTION, n_interface, cell_center, cell_current, dist,
      &dist_scaled, &dot_val, &dot_val_scaled);
  // double bot = sq(dot_val_scaled);
  // Malan Eq (27)
  // double bot = sq(dist_scaled);
  // double bot = 1.0;
  double bot = dist_scaled;

  double eps_col = 0.0;

  if (bot > 1.0e-10) {
    eps_col = fabs(collin) / bot;
  }
  return eps_col;
}

void update_interfacial_point(scalar interfacial_point) {
  foreach () {
    if (interfacial(point, f)) {
      interfacial_point[] = 1.0;
    } else {
      interfacial_point[] = 0.0;
    }
  }
  return;
}
void update_n_vec_interfacial_point(vector n_vec, scalar interfacial_point) {
  compute_normal(f, n_vec);
  update_interfacial_point(interfacial_point);
  return;
}

void clamp_step() {
  foreach () {
    f[] = (f[] < F_ERR ? 0. : f[] > 1. - F_ERR ? 1. : f[]);
  }
  return;
}
