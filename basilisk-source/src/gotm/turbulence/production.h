// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/production.F90

extern void production_ (
  integer * nlev,
  realtype * nn,
  realtype * ss,
  realtype * xp
);
static inline void turbulence_production (
  integer * nlev,
  realtype * nn,
  realtype * ss,
  realtype * xp) {
  production_ (nlev, nn, ss, xp);
}
