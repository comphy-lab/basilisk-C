// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/tkeeq.F90

extern void tkeeq_ (
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
static inline void turbulence_tkeeq (
  integer * nlev,
  realtype * dt,
  realtype * u_taus,
  realtype * u_taub,
  realtype * z0s,
  realtype * z0b,
  realtype * h,
  realtype * nn,
  realtype * ss) {
  tkeeq_ (nlev, dt, u_taus, u_taub, z0s, z0b, h, nn, ss);
}
