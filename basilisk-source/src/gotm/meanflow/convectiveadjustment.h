// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/convectiveadjustment.F90

extern void convectiveadjustment_ (
  integer * nlev,
  realtype * num,
  realtype * nuh,
  realtype * const_num,
  realtype * const_nuh,
  integer * buoy_method,
  realtype * g,
  realtype * rho_0
);
static inline void meanflow_convectiveadjustment (
  integer * nlev,
  realtype * num,
  realtype * nuh,
  realtype * const_num,
  realtype * const_nuh,
  integer * buoy_method,
  realtype * g,
  realtype * rho_0) {
  convectiveadjustment_ (nlev, num, nuh, const_num, const_nuh, buoy_method, g, rho_0);
}
