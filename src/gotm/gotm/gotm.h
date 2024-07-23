// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/gotm/gotm.F90

extern void __gotm_MOD_init_gotm (
  char * t1,
  char * t2,
  long int _slt1,
  long int _slt2
);
static inline void gotm_init_gotm (
  char * t1,
  char * t2) {
  __gotm_MOD_init_gotm (t1, t2, strlencheck(t1), strlencheck(t2));
}

extern void __gotm_MOD_time_loop (

);
static inline void gotm_time_loop (void) {
  __gotm_MOD_time_loop();
}

extern void __gotm_MOD_clean_up (

);
static inline void gotm_clean_up (void) {
  __gotm_MOD_clean_up();
}

extern void __gotm_MOD_print_state (

);
static inline void gotm_print_state (void) {
  __gotm_MOD_print_state();
}

extern void __gotm_MOD_setup_restart (

);
static inline void gotm_setup_restart (void) {
  __gotm_MOD_setup_restart();
}

extern void __gotm_MOD_read_restart (
  logical * restart_allow_missing_variable
);
static inline void gotm_read_restart (
  logical * restart_allow_missing_variable) {
  __gotm_MOD_read_restart (restart_allow_missing_variable);
}
