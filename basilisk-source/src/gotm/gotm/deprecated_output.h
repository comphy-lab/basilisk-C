// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/gotm/deprecated_output.F90

extern void deprecated_output_ (
  integer * namlst,
  char * title,
  realtype * dt,
  logical * list_fields,
  long int _sltitle
);
static inline void gotm_deprecated_output (
  integer * namlst,
  char * title,
  realtype * dt,
  logical * list_fields) {
  deprecated_output_ (namlst, title, dt, list_fields, strlencheck(title));
}
