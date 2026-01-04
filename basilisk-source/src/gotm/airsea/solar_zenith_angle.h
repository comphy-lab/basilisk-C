// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/solar_zenith_angle.F90

extern realtype solar_zenith_angle_ (
  integer * yday,
  realtype * hh,
  realtype * dlon,
  realtype * dlat
);
static inline realtype airsea_solar_zenith_angle (
  integer * yday,
  realtype * hh,
  realtype * dlon,
  realtype * dlat) {
  return solar_zenith_angle_ (yday, hh, dlon, dlat);
}
