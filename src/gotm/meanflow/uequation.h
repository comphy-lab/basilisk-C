// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/uequation.F90

extern void uequation_ (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * tx,
  realtype * num,
  realtype * gamu,
  integer * method
);
static inline void meanflow_uequation (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * tx,
  realtype * num,
  realtype * gamu,
  integer * method) {
  uequation_ (nlev, dt, cnpar, tx, num, gamu, method);
}
