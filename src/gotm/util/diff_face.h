// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/diff_face.F90

extern void diff_face_ (
  integer * n,
  realtype * dt,
  realtype * cnpar,
  realtype * h,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  realtype * nuy,
  realtype * lsour,
  realtype * qsour,
  realtype * y
);
static inline void util_diff_face (
  integer * n,
  realtype * dt,
  realtype * cnpar,
  realtype * h,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  realtype * nuy,
  realtype * lsour,
  realtype * qsour,
  realtype * y) {
  diff_face_ (n, dt, cnpar, h, bcup, bcdw, yup, ydw, nuy, lsour, qsour, y);
}
