
#ifndef SCALE_LS
#define SCALE_LS (1.0)
#endif

double sign2(double x) { return (x > 0. ? 1. : x < 0 ? -1. : 0.); }
double my_minmod(double a, double b) {
  if (a * b > 0) {
    if (fabs(a) < fabs(b)) return a;
    if (fabs(a) > fabs(b)) return b;
  }
  return 0;
}

foreach_dimension()
static inline double WENOdiff_x(Point point, scalar s, int i) {
  double s1 = (s[2 * i, 0, 0] + s[] - 2 * s[i, 0, 0]) / (Delta*SCALE_LS);
  double s2 = (s[1, 0, 0] + s[-1, 0, 0] - 2 * s[]) / (Delta*SCALE_LS);
  return i * ((s[i, 0, 0] - s[]) / (Delta*SCALE_LS) - my_minmod(s1, s2) / 2.);
}

void prehamil(Point point, coord* grapl, coord* gramin, scalar s) {
  foreach_dimension() {
    grapl->x = WENOdiff_x(point, s, 1);
    gramin->x = WENOdiff_x(point, s, -1);
  }
}

double hamiltonian(Point point, scalar s0, coord grapl, coord gramin) {
  double hamil = 0;
  if (s0[] > 0) {
    foreach_dimension() {
      double a = min(0., grapl.x);
      double b = max(0., gramin.x);
      hamil += max(sq(a), sq(b));
    }
    return sqrt(hamil);
  } else {
    foreach_dimension() {
      double a = max(0., grapl.x);
      double b = min(0., gramin.x);
      hamil += max(sq(a), sq(b));
    }
  }
  return sqrt(hamil);
}

foreach_dimension()
static inline double root_x(Point point, scalar s, double eps, int dir) {
  // dir == 1 or -1 offsets the position of the interface
  double phixx =
    my_minmod(s[2 * dir] + s[] - 2 * s[dir], s[1] + s[-1] - 2 * s[]);
  if (fabs(phixx) > eps) {
    double D = sq(phixx / 2. - s[] - s[dir]) - 4 * s[] * s[dir];
    // fprintf(stderr, "%g %g %g\n", D, phixx, eps);
    return 1 / 2. +
           (s[] - s[dir] - sign2(s[] - s[dir]) * sqrt(D)) / phixx;
  } else {
    double bot = (s[] - s[dir]);
    if (fabs(bot) < 1.0e-10) {
      bot = 1.0e-10;
    }
    // return s[] / (s[] - s[dir]);
    return s[] / bot;
  }
}

double ForwardEuler(scalar dist, scalar temp, scalar dist0, double dt) {
  double res = 0.;

  foreach (reduction(max:res)) {
    double delt = 0.;
    double flag = 1.;
    foreach_dimension() {
      flag = min(flag, dist0[-1] * dist0[]);
      flag = min(flag, dist0[1] * dist0[]);
    }

    coord grapl, gramin;
    prehamil(point, &grapl, &gramin, temp);

    if (flag < 0.) {  // the cell contains the interface
      double size = 1.e10;
      foreach_dimension() {
        if (dist0[] * dist0[1] < 0) {
          double dx = (Delta*SCALE_LS) * root_x(point, dist0, 1.e-10, 1);
          double sxx1 = (temp[2] + temp[] - 2 * temp[1]) / sq((Delta*SCALE_LS));
          double sxx2 = (temp[1] + temp[-1] - 2 * temp[]) / sq((Delta*SCALE_LS));
          if (dx != 0.)
            grapl.x = -temp[] / dx - dx * my_minmod(sxx1, sxx2) / 2.;
          else
            grapl.x = 0.;
          size = min(size, dx);
        }
        if (dist0[] * dist0[-1] < 0) {
            double dx = (Delta*SCALE_LS) * root_x(point, dist0, 1.e-10, -1);
            double sxx2 = (temp[1] + temp[-1] - 2 * temp[]) / sq((Delta*SCALE_LS));
            double sxx3 =
              (temp[-2] + temp[0] - 2 * temp[-1]) / sq((Delta*SCALE_LS));
            if (dx != 0.)
              gramin.x = temp[] / dx + dx * my_minmod(sxx3, sxx2) / 2.;
            else
              gramin.x = 0.;
            size = min(size, dx);
        }
      }
      delt = sign2(dist0[]) * min(dt, fabs(size) / 2.) *
             (hamiltonian(point, dist0, grapl, gramin) - 1);
      dist[] -= delt;
    } else {
      delt = sign2(dist0[]) * (hamiltonian(point, dist0, grapl, gramin) - 1);
      dist[] -= dt * delt;
    }
    res = max(res, fabs(delt));
  }

  restriction({dist});

  return res;
}

struct LS_reinit {
  scalar dist;
  double dt;
  int it_max;
};

int LS_reinit(struct LS_reinit p) {
  scalar dist = p.dist;

  double dt = p.dt;       // standard timestep (0.5*Delta)
  int it_max = p.it_max;  // maximum number of iteration (100)

  if (dt == 0) dt = 0.5 * (L0 / (1 << grid->maxdepth))*SCALE_LS;

  if (it_max == 0) it_max = 5;

  vector gr_LS[];
  int i;

  double eps = dt * 1.e-6;

  scalar dist0[];
  foreach () { dist0[] = dist[]; }

  for (i = 1; i <= it_max; i++) {
    double res = 0;

    scalar temp[], temp1[], temp2[];
    foreach () {
      temp[] = dist[];
      temp1[] = dist[];
    }
    //
    ForwardEuler(temp1, temp, dist0, dt);
    //
    foreach () { temp2[] = temp1[]; }
    //
    ForwardEuler(temp2, temp1, dist0, dt);
    //
    foreach () {
      temp1[] = 3. / 4 * dist[] + temp2[] / 4.;
      temp2[] = temp1[];
    }
    //

    ForwardEuler(temp2, temp1, dist0, dt);

    foreach (reduction(max:res)) {
      res = max(res, 2. / 3. * fabs(dist[] - temp2[]));
      dist[] = dist[] / 3. + temp2[] * 2. / 3.;
    }
    restriction({dist});
    if (res < eps) {
      return i;
    }
  }
  return it_max;
}
