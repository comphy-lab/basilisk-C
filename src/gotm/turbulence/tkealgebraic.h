// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/tkealgebraic.F90

extern void tkealgebraic_ (
  integer * nlev,
  realtype * u_taus,
  realtype * u_taub,
  realtype * nn,
  realtype * ss
);
static inline void turbulence_tkealgebraic (
  integer * nlev,
  realtype * u_taus,
  realtype * u_taub,
  realtype * nn,
  realtype * ss) {
  tkealgebraic_ (nlev, u_taus, u_taub, nn, ss);
}
