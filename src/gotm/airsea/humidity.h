// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/humidity.F90

extern void humidity_ (
  integer * hum_method,
  realtype * hum,
  realtype * airp,
  realtype * tw,
  realtype * ta
);
static inline void airsea_humidity (
  integer * hum_method,
  realtype * hum,
  realtype * airp,
  realtype * tw,
  realtype * ta) {
  humidity_ (hum_method, hum, airp, tw, ta);
}
