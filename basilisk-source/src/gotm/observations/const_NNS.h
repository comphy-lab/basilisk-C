// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/observations/const_NNS.F90

extern void const_nns_ (
  integer * nlev,
  realtype * z,
  realtype * s_top,
  realtype * t_const,
  realtype * nn,
  realtype * gravity,
  realtype * rho_0,
  realtype * s
);
static inline void observations_const_nns (
  integer * nlev,
  realtype * z,
  realtype * s_top,
  realtype * t_const,
  realtype * nn,
  realtype * gravity,
  realtype * rho_0,
  realtype * s) {
  const_nns_ (nlev, z, s_top, t_const, nn, gravity, rho_0, s);
}
