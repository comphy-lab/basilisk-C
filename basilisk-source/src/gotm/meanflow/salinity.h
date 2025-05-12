// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/salinity.F90

extern void salinity_ (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * nus,
  realtype * gams
);
static inline void meanflow_salinity (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * nus,
  realtype * gams) {
  salinity_ (nlev, dt, cnpar, nus, gams);
}
