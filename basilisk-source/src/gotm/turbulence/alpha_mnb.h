// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/alpha_mnb.F90

extern void alpha_mnb_ (
  integer * nlev,
  realtype * nn,
  realtype * ss
);
static inline void turbulence_alpha_mnb (
  integer * nlev,
  realtype * nn,
  realtype * ss) {
  alpha_mnb_ (nlev, nn, ss);
}
