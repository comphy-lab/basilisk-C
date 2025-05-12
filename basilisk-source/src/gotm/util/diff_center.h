// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/diff_center.F90

extern void diff_center_ (
  integer * n,
  realtype * dt,
  realtype * cnpar,
  integer * posconc,
  realtype * h,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  realtype * nuy,
  realtype * lsour,
  realtype * qsour,
  realtype * taur,
  realtype * yobs,
  realtype * y
);
static inline void util_diff_center (
  integer * n,
  realtype * dt,
  realtype * cnpar,
  integer * posconc,
  realtype * h,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  realtype * nuy,
  realtype * lsour,
  realtype * qsour,
  realtype * taur,
  realtype * yobs,
  realtype * y) {
  diff_center_ (n, dt, cnpar, posconc, h, bcup, bcdw, yup, ydw, nuy, lsour, qsour, taur, yobs, y);
}
