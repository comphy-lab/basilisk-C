// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/gotm/diagnostics.F90
extern realtype __diagnostics_MOD_ekin;
#define diagnostics_ekin __diagnostics_MOD_ekin
extern realtype __diagnostics_MOD_epot;
#define diagnostics_epot __diagnostics_MOD_epot
extern realtype __diagnostics_MOD_eturb;
#define diagnostics_eturb __diagnostics_MOD_eturb
extern realtype __diagnostics_MOD_taux;
#define diagnostics_taux __diagnostics_MOD_taux
extern realtype __diagnostics_MOD_tauy;
#define diagnostics_tauy __diagnostics_MOD_tauy
extern integer __diagnostics_MOD_mld_method;
#define diagnostics_mld_method __diagnostics_MOD_mld_method
extern realtype __diagnostics_MOD_mld_surf;
#define diagnostics_mld_surf __diagnostics_MOD_mld_surf
extern realtype __diagnostics_MOD_mld_bott;
#define diagnostics_mld_bott __diagnostics_MOD_mld_bott

extern void __diagnostics_MOD_init_diagnostics (
  integer * nlev
);
static inline void diagnostics_init_diagnostics (
  integer * nlev) {
  __diagnostics_MOD_init_diagnostics (nlev);
}

extern void __diagnostics_MOD_do_diagnostics (
  integer * nlev
);
static inline void diagnostics_do_diagnostics (
  integer * nlev) {
  __diagnostics_MOD_do_diagnostics (nlev);
}

extern void __diagnostics_MOD_clean_diagnostics (

);
static inline void diagnostics_clean_diagnostics (void) {
  __diagnostics_MOD_clean_diagnostics();
}
realtype diagnostics_get_global (const char * name) {
  if (!strcmp (name, "ekin"))
    return diagnostics_ekin;
  if (!strcmp (name, "epot"))
    return diagnostics_epot;
  if (!strcmp (name, "eturb"))
    return diagnostics_eturb;
  if (!strcmp (name, "taux"))
    return diagnostics_taux;
  if (!strcmp (name, "tauy"))
    return diagnostics_tauy;
  if (!strcmp (name, "mld_method"))
    return diagnostics_mld_method;
  if (!strcmp (name, "mld_surf"))
    return diagnostics_mld_surf;
  if (!strcmp (name, "mld_bott"))
    return diagnostics_mld_bott;
  return HUGE;
}
