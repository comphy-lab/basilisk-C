// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/shear.F90

extern void shear_ (
  integer * nlev,
  realtype * cnpar
);
static inline void meanflow_shear (
  integer * nlev,
  realtype * cnpar) {
  shear_ (nlev, cnpar);
}
