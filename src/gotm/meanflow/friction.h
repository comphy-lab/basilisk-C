// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/meanflow/friction.F90

extern void friction_ (
  realtype * kappa,
  realtype * avmolu,
  realtype * tx,
  realtype * ty
);
static inline void meanflow_friction (
  realtype * kappa,
  realtype * avmolu,
  realtype * tx,
  realtype * ty) {
  friction_ (kappa, avmolu, tx, ty);
}
