// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/coriolis.F90

extern void coriolis_ (
  integer * nlev,
  realtype * dt
);
static inline void meanflow_coriolis (
  integer * nlev,
  realtype * dt) {
  coriolis_ (nlev, dt);
}
