// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/ode_solvers.F90

extern void findp_bisection_ (
  integer * numc,
  realtype * cc,
  realtype * derivative,
  realtype * dt,
  realtype * accuracy,
  realtype * pi
);
static inline void util_findp_bisection (
  integer * numc,
  realtype * cc,
  realtype * derivative,
  realtype * dt,
  realtype * accuracy,
  realtype * pi) {
  findp_bisection_ (numc, cc, derivative, dt, accuracy, pi);
}

extern void matrix_ (
  integer * n,
  realtype * a,
  realtype * r,
  realtype * c
);
static inline void util_matrix (
  integer * n,
  realtype * a,
  realtype * r,
  realtype * c) {
  matrix_ (n, a, r, c);
}
