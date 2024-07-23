// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/lengthscaleeq.F90

extern void lengthscaleeq_ (
  integer * nlev,
  realtype * dt,
  realtype * depth,
  realtype * u_taus,
  realtype * u_taub,
  realtype * z0s,
  realtype * z0b,
  realtype * h,
  realtype * nn,
  realtype * ss
);
static inline void turbulence_lengthscaleeq (
  integer * nlev,
  realtype * dt,
  realtype * depth,
  realtype * u_taus,
  realtype * u_taub,
  realtype * z0s,
  realtype * z0b,
  realtype * h,
  realtype * nn,
  realtype * ss) {
  lengthscaleeq_ (nlev, dt, depth, u_taus, u_taub, z0s, z0b, h, nn, ss);
}
