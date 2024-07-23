// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/adv_center.F90

extern void adv_center_ (
  integer * n,
  realtype * dt,
  realtype * h,
  realtype * ho,
  realtype * ww,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  integer * method,
  integer * mode,
  realtype * y
);
static inline void util_adv_center (
  integer * n,
  realtype * dt,
  realtype * h,
  realtype * ho,
  realtype * ww,
  integer * bcup,
  integer * bcdw,
  realtype * yup,
  realtype * ydw,
  integer * method,
  integer * mode,
  realtype * y) {
  adv_center_ (n, dt, h, ho, ww, bcup, bcdw, yup, ydw, method, mode, y);
}
