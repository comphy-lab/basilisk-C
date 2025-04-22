// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/variances.F90

extern void variances_ (
  integer * nlev,
  realtype * ssu,
  realtype * ssv
);
static inline void turbulence_variances (
  integer * nlev,
  realtype * ssu,
  realtype * ssv) {
  variances_ (nlev, ssu, ssv);
}
