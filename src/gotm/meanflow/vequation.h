// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/vequation.F90

extern void vequation_ (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * ty,
  realtype * num,
  realtype * gamv,
  integer * method
);
static inline void meanflow_vequation (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * ty,
  realtype * num,
  realtype * gamv,
  integer * method) {
  vequation_ (nlev, dt, cnpar, ty, num, gamv, method);
}
