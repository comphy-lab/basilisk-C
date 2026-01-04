// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/algebraiclength.F90

extern void algebraiclength_ (
  integer * method,
  integer * nlev,
  realtype * z0b,
  realtype * z0s,
  realtype * depth,
  realtype * h,
  realtype * nn
);
static inline void turbulence_algebraiclength (
  integer * method,
  integer * nlev,
  realtype * z0b,
  realtype * z0s,
  realtype * depth,
  realtype * h,
  realtype * nn) {
  algebraiclength_ (method, nlev, z0b, z0s, depth, h, nn);
}
