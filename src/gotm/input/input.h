// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/input/input.F90

extern void __input_MOD_init_input (
  integer * n
);
static inline void input_init_input (
  integer * n) {
  __input_MOD_init_input (n);
}

extern void __input_MOD_register_input_1d (
  char * path,
  integer * icolumn,
  realtype_1d * data,
  char * name,
  realtype * scale_factor,
  long int _slpath,
  long int _slname
);
static inline void input_register_input_1d (
  char * path,
  integer * icolumn,
  realtype_1d * data,
  char * name,
  realtype * scale_factor) {
  __input_MOD_register_input_1d (path, icolumn, data, name, scale_factor, strlencheck(path), strlencheck(name));
}

extern void __input_MOD_register_input_0d (
  char * path,
  integer * icolumn,
  realtype * data,
  char * name,
  realtype * scale_factor,
  realtype * add_offset,
  long int _slpath,
  long int _slname
);
static inline void input_register_input_0d (
  char * path,
  integer * icolumn,
  realtype * data,
  char * name,
  realtype * scale_factor,
  realtype * add_offset) {
  __input_MOD_register_input_0d (path, icolumn, data, name, scale_factor, add_offset, strlencheck(path), strlencheck(name));
}

extern void __input_MOD_do_input (
  integer * jul,
  integer * secs,
  integer * nlev,
  realtype_1d * z
);
static inline void input_do_input (
  integer * jul,
  integer * secs,
  integer * nlev,
  realtype_1d * z) {
  __input_MOD_do_input (jul, secs, nlev, z);
}

extern void __input_MOD_initialize_profile_file (
  type_profile_file * info,
  integer * nlev
);
static inline void input_initialize_profile_file (
  type_profile_file * info,
  integer * nlev) {
  __input_MOD_initialize_profile_file (info, nlev);
}

extern void __input_MOD_get_observed_profiles (
  type_profile_file * info,
  integer * jul,
  integer * secs,
  integer * nlev,
  realtype * z
);
static inline void input_get_observed_profiles (
  type_profile_file * info,
  integer * jul,
  integer * secs,
  integer * nlev,
  realtype * z) {
  __input_MOD_get_observed_profiles (info, jul, secs, nlev, z);
}

extern void __input_MOD_initialize_timeseries_file (
  type_timeseries_file * info
);
static inline void input_initialize_timeseries_file (
  type_timeseries_file * info) {
  __input_MOD_initialize_timeseries_file (info);
}

extern void __input_MOD_get_observed_scalars (
  type_timeseries_file * info,
  integer * jul,
  integer * secs
);
static inline void input_get_observed_scalars (
  type_timeseries_file * info,
  integer * jul,
  integer * secs) {
  __input_MOD_get_observed_scalars (info, jul, secs);
}

extern void __input_MOD_close_input (

);
static inline void input_close_input (void) {
  __input_MOD_close_input();
}

extern void __input_MOD_read_obs (
  integer * unit,
  integer * yy,
  integer * mm,
  integer * dd,
  integer * hh,
  integer * min,
  integer * ss,
  integer * n,
  realtype_1d * obs,
  integer * ios,
  integer * line
);
static inline void input_read_obs (
  integer * unit,
  integer * yy,
  integer * mm,
  integer * dd,
  integer * hh,
  integer * min,
  integer * ss,
  integer * n,
  realtype_1d * obs,
  integer * ios,
  integer * line) {
  __input_MOD_read_obs (unit, yy, mm, dd, hh, min, ss, n, obs, ios, line);
}

extern void __input_MOD_read_profiles (
  integer * unit,
  integer * nlev,
  integer * cols,
  integer * yy,
  integer * mm,
  integer * dd,
  integer * hh,
  integer * min,
  integer * ss,
  realtype_1d * z,
  realtype * profiles,
  integer * lines,
  integer * ios
);
static inline void input_read_profiles (
  integer * unit,
  integer * nlev,
  integer * cols,
  integer * yy,
  integer * mm,
  integer * dd,
  integer * hh,
  integer * min,
  integer * ss,
  realtype_1d * z,
  realtype * profiles,
  integer * lines,
  integer * ios) {
  __input_MOD_read_profiles (unit, nlev, cols, yy, mm, dd, hh, min, ss, z, profiles, lines, ios);
}

extern void __input_MOD_fatal_error (
  char * location,
  char * error,
  long int _sllocation,
  long int _slerror
);
static inline void input_fatal_error (
  char * location,
  char * error) {
  __input_MOD_fatal_error (location, error, strlencheck(location), strlencheck(error));
}
