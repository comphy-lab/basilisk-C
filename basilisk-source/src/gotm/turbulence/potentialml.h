// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/potentialml.F90

extern void potentialml_ (
  integer * nlev,
  realtype * z0b,
  realtype * z0s,
  realtype * h,
  realtype * depth,
  realtype * nn
);
static inline void turbulence_potentialml (
  integer * nlev,
  realtype * z0b,
  realtype * z0s,
  realtype * h,
  realtype * depth,
  realtype * nn) {
  potentialml_ (nlev, z0b, z0s, h, depth, nn);
}
