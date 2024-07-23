// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/eqstate.F90
extern integer __eqstate_MOD_eq_state_method;
#define eqstate_eq_state_method __eqstate_MOD_eq_state_method
extern integer __eqstate_MOD_eq_state_mode;
#define eqstate_eq_state_mode __eqstate_MOD_eq_state_mode
extern realtype __eqstate_MOD_t0;
#define eqstate_t0 __eqstate_MOD_t0
extern realtype __eqstate_MOD_s0;
#define eqstate_s0 __eqstate_MOD_s0
extern realtype __eqstate_MOD_p0;
#define eqstate_p0 __eqstate_MOD_p0
extern realtype __eqstate_MOD_dtr0;
#define eqstate_dtr0 __eqstate_MOD_dtr0
extern realtype __eqstate_MOD_dsr0;
#define eqstate_dsr0 __eqstate_MOD_dsr0

extern void __eqstate_MOD_init_eqstate (
  integer * namlst
);
static inline void eqstate_init_eqstate (
  integer * namlst) {
  __eqstate_MOD_init_eqstate (namlst);
}

extern realtype __eqstate_MOD_eqstate1 (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0
);
static inline realtype eqstate_eqstate1 (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0) {
  return __eqstate_MOD_eqstate1 (s, t, p, g, rho_0);
}

extern realtype __eqstate_MOD_eos_alpha (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0
);
static inline realtype eqstate_eos_alpha (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0) {
  return __eqstate_MOD_eos_alpha (s, t, p, g, rho_0);
}

extern realtype __eqstate_MOD_eos_beta (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0
);
static inline realtype eqstate_eos_beta (
  realtype * s,
  realtype * t,
  realtype * p,
  realtype * g,
  realtype * rho_0) {
  return __eqstate_MOD_eos_beta (s, t, p, g, rho_0);
}

extern realtype __eqstate_MOD_unesco (
  realtype * s,
  realtype * t,
  realtype * p,
  logical * unpress
);
static inline realtype eqstate_unesco (
  realtype * s,
  realtype * t,
  realtype * p,
  logical * unpress) {
  return __eqstate_MOD_unesco (s, t, p, unpress);
}

extern realtype __eqstate_MOD_rho_feistel (
  realtype * s,
  realtype * th,
  realtype * p,
  logical * unpress
);
static inline realtype eqstate_rho_feistel (
  realtype * s,
  realtype * th,
  realtype * p,
  logical * unpress) {
  return __eqstate_MOD_rho_feistel (s, th, p, unpress);
}
realtype eqstate_get_global (const char * name) {
  if (!strcmp (name, "eq_state_method"))
    return eqstate_eq_state_method;
  if (!strcmp (name, "eq_state_mode"))
    return eqstate_eq_state_mode;
  if (!strcmp (name, "t0"))
    return eqstate_t0;
  if (!strcmp (name, "s0"))
    return eqstate_s0;
  if (!strcmp (name, "p0"))
    return eqstate_p0;
  if (!strcmp (name, "dtr0"))
    return eqstate_dtr0;
  if (!strcmp (name, "dsr0"))
    return eqstate_dsr0;
  return HUGE;
}
