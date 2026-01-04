// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/albedo_water.F90

extern realtype albedo_water_ (
  integer * method,
  realtype * zenith_angle,
  integer * yday
);
static inline realtype airsea_albedo_water (
  integer * method,
  realtype * zenith_angle,
  integer * yday) {
  return albedo_water_ (method, zenith_angle, yday);
}

extern realtype albedo_payne_ (
  realtype * zen
);
static inline realtype airsea_albedo_payne (
  realtype * zen) {
  return albedo_payne_ (zen);
}

extern realtype albedo_cogley_ (
  realtype * zen,
  integer * yd
);
static inline realtype airsea_albedo_cogley (
  realtype * zen,
  integer * yd) {
  return albedo_cogley_ (zen, yd);
}
