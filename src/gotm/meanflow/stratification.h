// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/stratification.F90

extern void stratification_ (
  integer * nlev,
  integer * buoy_method,
  realtype * dt,
  realtype * cnpar,
  realtype * nub,
  realtype * gamb
);
static inline void meanflow_stratification (
  integer * nlev,
  integer * buoy_method,
  realtype * dt,
  realtype * cnpar,
  realtype * nub,
  realtype * gamb) {
  stratification_ (nlev, buoy_method, dt, cnpar, nub, gamb);
}
