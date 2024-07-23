// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/wequation.F90

extern void wequation_ (
  integer * nlev,
  realtype * dt
);
static inline void meanflow_wequation (
  integer * nlev,
  realtype * dt) {
  wequation_ (nlev, dt);
}
