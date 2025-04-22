// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/temperature.F90

extern void temperature_ (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * i_0,
  realtype * heat,
  realtype * nuh,
  realtype * gamh,
  realtype * rad
);
static inline void meanflow_temperature (
  integer * nlev,
  realtype * dt,
  realtype * cnpar,
  realtype * i_0,
  realtype * heat,
  realtype * nuh,
  realtype * gamh,
  realtype * rad) {
  temperature_ (nlev, dt, cnpar, i_0, heat, nuh, gamh, rad);
}
