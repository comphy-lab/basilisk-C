// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/tridiagonal.F90

extern void __mtridiagonal_MOD_init_tridiagonal (
  integer * n
);
static inline void mtridiagonal_init_tridiagonal (
  integer * n) {
  __mtridiagonal_MOD_init_tridiagonal (n);
}

extern void __mtridiagonal_MOD_tridiagonal (
  integer * n,
  integer * fi,
  integer * lt,
  realtype * value
);
static inline void mtridiagonal_tridiagonal (
  integer * n,
  integer * fi,
  integer * lt,
  realtype * value) {
  __mtridiagonal_MOD_tridiagonal (n, fi, lt, value);
}

extern void __mtridiagonal_MOD_clean_tridiagonal (

);
static inline void mtridiagonal_clean_tridiagonal (void) {
  __mtridiagonal_MOD_clean_tridiagonal();
}
