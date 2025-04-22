// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/buoyancy.F90

extern void buoyancy_ (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * nub,
  realtype * gamb
);
static inline void meanflow_buoyancy (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * nub,
  realtype * gamb) {
  buoyancy_ (nlev, dt, cnpar, nub, gamb);
}
