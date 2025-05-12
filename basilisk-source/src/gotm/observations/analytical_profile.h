// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/observations/analytical_profile.F90

extern void analytical_profile_ (
  integer * nlev,
  realtype * z,
  realtype * z1,
  realtype * v1,
  realtype * z2,
  realtype * v2,
  realtype * prof
);
static inline void observations_analytical_profile (
  integer * nlev,
  realtype * z,
  realtype * z1,
  realtype * v1,
  realtype * z2,
  realtype * v2,
  realtype * prof) {
  analytical_profile_ (nlev, z, z1, v1, z2, v2, prof);
}
