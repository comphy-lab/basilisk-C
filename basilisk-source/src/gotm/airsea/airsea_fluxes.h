// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/airsea_fluxes.F90

extern void airsea_fluxes_ (
  integer * method,
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
static inline void airsea_airsea_fluxes (
  integer * method,
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
  airsea_fluxes_ (method, sst, airt, u10, v10, precip, evap, taux, tauy, qe, qh);
}
