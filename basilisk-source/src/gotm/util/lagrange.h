// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/lagrange.F90

extern void lagrange_ (
  realtype * nlev,
  realtype * dt,
  realtype * zlev,
  realtype * nuh,
  realtype * w,
  realtype * npar,
  logical * active,
  integer * zi,
  realtype * zp
);
static inline void util_lagrange (
  realtype * nlev,
  realtype * dt,
  realtype * zlev,
  realtype * nuh,
  realtype * w,
  realtype * npar,
  logical * active,
  integer * zi,
  realtype * zp) {
  lagrange_ (nlev, dt, zlev, nuh, w, npar, active, zi, zp);
}
