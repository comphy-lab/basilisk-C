// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/observations/const_NNT.F90

extern void const_nnt_ (
  integer * nlev,
  realtype * z,
  realtype * t_top,
  realtype * s_const,
  realtype * nn,
  realtype * gravity,
  realtype * rho_0,
  realtype * t
);
static inline void observations_const_nnt (
  integer * nlev,
  realtype * z,
  realtype * t_top,
  realtype * s_const,
  realtype * nn,
  realtype * gravity,
  realtype * rho_0,
  realtype * t) {
  const_nnt_ (nlev, z, t_top, s_const, nn, gravity, rho_0, t);
}
