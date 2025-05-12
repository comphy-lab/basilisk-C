// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/ispralength.F90

extern void ispralength_ (
  integer * nlev,
  realtype * nn,
  realtype * h,
  realtype * depth
);
static inline void turbulence_ispralength (
  integer * nlev,
  realtype * nn,
  realtype * h,
  realtype * depth) {
  ispralength_ (nlev, nn, h, depth);
}
