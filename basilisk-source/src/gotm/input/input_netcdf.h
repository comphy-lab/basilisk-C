// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/input/input_netcdf.F90

extern void __input_netcdf_MOD_init_input_netcdf (
  integer * n
);
static inline void input_netcdf_init_input_netcdf (
  integer * n) {
  __input_netcdf_MOD_init_input_netcdf (n);
}

extern void __input_netcdf_MOD_read_restart_data (
  char * var_name,
  logical * allow_missing_variable,
  realtype * data_0d,
  realtype_1d * data_1d,
  long int _slvar_name
);
static inline void input_netcdf_read_restart_data (
  char * var_name,
  logical * allow_missing_variable,
  realtype * data_0d,
  realtype_1d * data_1d) {
  __input_netcdf_MOD_read_restart_data (var_name, allow_missing_variable, data_0d, data_1d, strlencheck(var_name));
}

extern void __input_netcdf_MOD_handle_err (
  integer * ierr,
  char * msg,
  long int _slmsg
);
static inline void input_netcdf_handle_err (
  integer * ierr,
  char * msg) {
  __input_netcdf_MOD_handle_err (ierr, msg, strlencheck(msg));
}
