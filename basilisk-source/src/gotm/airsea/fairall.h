// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/fairall.F90

extern void fairall_ (
  realtype * sst,
  realtype * airt,
  realtype * u10,
  realtype * v10,
  realtype * precip,
  realtype * evap,
  realtype * taux,
  realtype * tauy,
  realtype * qe,
  realtype * qh
);
static inline void airsea_fairall (
  realtype * sst,
  realtype * airt,
  realtype * u10,
  realtype * v10,
  realtype * precip,
  realtype * evap,
  realtype * taux,
  realtype * tauy,
  realtype * qe,
  realtype * qh) {
  fairall_ (sst, airt, u10, v10, precip, evap, taux, tauy, qe, qh);
}

extern realtype psi_ (
  integer * iflag,
  realtype * zol
);
static inline realtype airsea_psi (
  integer * iflag,
  realtype * zol) {
  return psi_ (iflag, zol);
}
