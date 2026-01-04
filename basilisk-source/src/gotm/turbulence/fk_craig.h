// Generated automatically by gotm.awk from
// /home/popinet/local/src/GOTM-5.2.1/src/turbulence/fk_craig.F90

extern realtype fk_craig_ (
  realtype * u_tau,
  realtype * eta
);
static inline realtype turbulence_fk_craig (
  realtype * u_tau,
  realtype * eta) {
  return fk_craig_ (u_tau, eta);
}
