// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/q2over2eq.F90

extern void q2over2eq_ (
  integer * nlev,
  realtype * dt,
  realtype * u_taus,
  realtype * u_taub,
  realtype * z0s,
  realtype * z0b,
  realtype * h,
  realtype * nn,
  realtype * ss
);
static inline void turbulence_q2over2eq (
  integer * nlev,
  realtype * dt,
  realtype * u_taus,
  realtype * u_taub,
  realtype * z0s,
  realtype * z0b,
  realtype * h,
  realtype * nn,
  realtype * ss) {
  q2over2eq_ (nlev, dt, u_taus, u_taub, z0s, z0b, h, nn, ss);
}
