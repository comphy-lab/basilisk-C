// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/internal_wave.F90

extern void internal_wave_ (
  integer * nlev,
  realtype * nn,
  realtype * ss
);
static inline void turbulence_internal_wave (
  integer * nlev,
  realtype * nn,
  realtype * ss) {
  internal_wave_ (nlev, nn, ss);
}
