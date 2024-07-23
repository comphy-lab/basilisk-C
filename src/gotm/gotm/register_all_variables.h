// Generated automatically by gotm.awk from
// /home/popinet/local/src/GOTM-5.2.1/src/gotm/register_all_variables.F90
extern type_field_manager __register_all_variables_MOD_fm;
#define register_all_variables_fm __register_all_variables_MOD_fm

extern void __register_all_variables_MOD_do_register_all_variables (
  realtype * lat,
  realtype * lon,
  integer * nlev
);
static inline void register_all_variables_do_register_all_variables (
  realtype * lat,
  realtype * lon,
  integer * nlev) {
  __register_all_variables_MOD_do_register_all_variables (lat, lon, nlev);
}

extern void __register_all_variables_MOD_register_coordinate_variables (
  realtype * lat,
  realtype * lon
);
static inline void register_all_variables_register_coordinate_variables (
  realtype * lat,
  realtype * lon) {
  __register_all_variables_MOD_register_coordinate_variables (lat, lon);
}

extern void __register_all_variables_MOD_register_airsea_variables (
  integer * nlev
);
static inline void register_all_variables_register_airsea_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_airsea_variables (nlev);
}

extern void __register_all_variables_MOD_register_observation_variables (
  integer * nlev
);
static inline void register_all_variables_register_observation_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_observation_variables (nlev);
}

extern void __register_all_variables_MOD_register_stream_variables (
  integer * nlev
);
static inline void register_all_variables_register_stream_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_stream_variables (nlev);
}

extern void __register_all_variables_MOD_register_meanflow_variables (
  integer * nlev
);
static inline void register_all_variables_register_meanflow_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_meanflow_variables (nlev);
}

extern void __register_all_variables_MOD_register_diagnostic_variables (
  integer * nlev
);
static inline void register_all_variables_register_diagnostic_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_diagnostic_variables (nlev);
}

extern void __register_all_variables_MOD_register_aaa_variables (
  integer * nlev
);
static inline void register_all_variables_register_aaa_variables (
  integer * nlev) {
  __register_all_variables_MOD_register_aaa_variables (nlev);
}
