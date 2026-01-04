// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/updategrid.F90

extern void updategrid_ (
  integer * nlev,
  realtype * dt,
  realtype * zeta
);
static inline void meanflow_updategrid (
  integer * nlev,
  realtype * dt,
  realtype * zeta) {
  updategrid_ (nlev, dt, zeta);
}
