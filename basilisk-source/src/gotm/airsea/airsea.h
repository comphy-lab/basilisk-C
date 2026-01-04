// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/airsea.F90
extern logical __airsea_MOD_calc_fluxes;
#define airsea_calc_fluxes __airsea_MOD_calc_fluxes
extern integer __airsea_MOD_hum_method;
#define airsea_hum_method __airsea_MOD_hum_method
extern realtype __airsea_MOD_u10;
#define airsea_u10 __airsea_MOD_u10
extern realtype __airsea_MOD_v10;
#define airsea_v10 __airsea_MOD_v10
extern realtype __airsea_MOD_airp;
#define airsea_airp __airsea_MOD_airp
extern realtype __airsea_MOD_airt;
#define airsea_airt __airsea_MOD_airt
extern realtype __airsea_MOD_hum;
#define airsea_hum __airsea_MOD_hum
extern realtype __airsea_MOD_cloud;
#define airsea_cloud __airsea_MOD_cloud
extern realtype __airsea_MOD_w;
#define airsea_w __airsea_MOD_w
extern realtype __airsea_MOD_i_0;
#define airsea_i_0 __airsea_MOD_i_0
extern realtype __airsea_MOD_albedo;
#define airsea_albedo __airsea_MOD_albedo
extern realtype __airsea_MOD_heat;
#define airsea_heat __airsea_MOD_heat
extern realtype __airsea_MOD_qe;
#define airsea_qe __airsea_MOD_qe
extern realtype __airsea_MOD_qh;
#define airsea_qh __airsea_MOD_qh
extern realtype __airsea_MOD_qb;
#define airsea_qb __airsea_MOD_qb
extern realtype __airsea_MOD_tx;
#define airsea_tx __airsea_MOD_tx
extern realtype __airsea_MOD_ty;
#define airsea_ty __airsea_MOD_ty
extern realtype __airsea_MOD_precip;
#define airsea_precip __airsea_MOD_precip
extern realtype __airsea_MOD_evap;
#define airsea_evap __airsea_MOD_evap
extern realtype __airsea_MOD_sst;
#define airsea_sst __airsea_MOD_sst
extern realtype __airsea_MOD_sst_obs;
#define airsea_sst_obs __airsea_MOD_sst_obs
extern realtype __airsea_MOD_sss;
#define airsea_sss __airsea_MOD_sss
extern realtype __airsea_MOD_ssu;
#define airsea_ssu __airsea_MOD_ssu
extern realtype __airsea_MOD_ssv;
#define airsea_ssv __airsea_MOD_ssv
extern realtype __airsea_MOD_int_precip;
#define airsea_int_precip __airsea_MOD_int_precip
extern realtype __airsea_MOD_int_evap;
#define airsea_int_evap __airsea_MOD_int_evap
extern realtype __airsea_MOD_int_fwf;
#define airsea_int_fwf __airsea_MOD_int_fwf
extern realtype __airsea_MOD_int_swr;
#define airsea_int_swr __airsea_MOD_int_swr
extern realtype __airsea_MOD_int_heat;
#define airsea_int_heat __airsea_MOD_int_heat
extern realtype __airsea_MOD_int_total;
#define airsea_int_total __airsea_MOD_int_total
extern realtype __airsea_MOD_bio_drag_scale;
#define airsea_bio_drag_scale __airsea_MOD_bio_drag_scale
extern realtype __airsea_MOD_bio_albedo;
#define airsea_bio_albedo __airsea_MOD_bio_albedo
extern integer __airsea_MOD_swr_method;
#define airsea_swr_method __airsea_MOD_swr_method
extern integer __airsea_MOD_albedo_method;
#define airsea_albedo_method __airsea_MOD_albedo_method
extern realtype __airsea_MOD_const_albedo;
#define airsea_const_albedo __airsea_MOD_const_albedo
extern integer __airsea_MOD_fluxes_method;
#define airsea_fluxes_method __airsea_MOD_fluxes_method
extern integer __airsea_MOD_back_radiation_method;
#define airsea_back_radiation_method __airsea_MOD_back_radiation_method
extern logical __airsea_MOD_rain_impact;
#define airsea_rain_impact __airsea_MOD_rain_impact
extern logical __airsea_MOD_calc_evaporation;
#define airsea_calc_evaporation __airsea_MOD_calc_evaporation

extern void __airsea_MOD_init_air_sea (
  integer * namlst,
  realtype * lat,
  realtype * lon
);
static inline void airsea_init_air_sea (
  integer * namlst,
  realtype * lat,
  realtype * lon) {
  __airsea_MOD_init_air_sea (namlst, lat, lon);
}

extern void __airsea_MOD_do_air_sea (
  integer * jul,
  integer * secs
);
static inline void airsea_do_air_sea (
  integer * jul,
  integer * secs) {
  __airsea_MOD_do_air_sea (jul, secs);
}

extern void __airsea_MOD_clean_air_sea (

);
static inline void airsea_clean_air_sea (void) {
  __airsea_MOD_clean_air_sea();
}

extern void __airsea_MOD_flux_from_meteo (
  integer * jul,
  integer * secs
);
static inline void airsea_flux_from_meteo (
  integer * jul,
  integer * secs) {
  __airsea_MOD_flux_from_meteo (jul, secs);
}

extern void __airsea_MOD_integrated_fluxes (
  realtype * dt
);
static inline void airsea_integrated_fluxes (
  realtype * dt) {
  __airsea_MOD_integrated_fluxes (dt);
}

extern void __airsea_MOD_set_sst (
  realtype * temp
);
static inline void airsea_set_sst (
  realtype * temp) {
  __airsea_MOD_set_sst (temp);
}

extern void __airsea_MOD_set_ssuv (
  realtype * uvel,
  realtype * vvel
);
static inline void airsea_set_ssuv (
  realtype * uvel,
  realtype * vvel) {
  __airsea_MOD_set_ssuv (uvel, vvel);
}

extern void __airsea_MOD_print_state_airsea (

);
static inline void airsea_print_state_airsea (void) {
  __airsea_MOD_print_state_airsea();
}
realtype airsea_get_global (const char * name) {
  if (!strcmp (name, "calc_fluxes"))
    return airsea_calc_fluxes;
  if (!strcmp (name, "hum_method"))
    return airsea_hum_method;
  if (!strcmp (name, "u10"))
    return airsea_u10;
  if (!strcmp (name, "v10"))
    return airsea_v10;
  if (!strcmp (name, "airp"))
    return airsea_airp;
  if (!strcmp (name, "airt"))
    return airsea_airt;
  if (!strcmp (name, "hum"))
    return airsea_hum;
  if (!strcmp (name, "cloud"))
    return airsea_cloud;
  if (!strcmp (name, "w"))
    return airsea_w;
  if (!strcmp (name, "i_0"))
    return airsea_i_0;
  if (!strcmp (name, "albedo"))
    return airsea_albedo;
  if (!strcmp (name, "heat"))
    return airsea_heat;
  if (!strcmp (name, "qe"))
    return airsea_qe;
  if (!strcmp (name, "qh"))
    return airsea_qh;
  if (!strcmp (name, "qb"))
    return airsea_qb;
  if (!strcmp (name, "tx"))
    return airsea_tx;
  if (!strcmp (name, "ty"))
    return airsea_ty;
  if (!strcmp (name, "precip"))
    return airsea_precip;
  if (!strcmp (name, "evap"))
    return airsea_evap;
  if (!strcmp (name, "sst"))
    return airsea_sst;
  if (!strcmp (name, "sst_obs"))
    return airsea_sst_obs;
  if (!strcmp (name, "sss"))
    return airsea_sss;
  if (!strcmp (name, "ssu"))
    return airsea_ssu;
  if (!strcmp (name, "ssv"))
    return airsea_ssv;
  if (!strcmp (name, "int_precip"))
    return airsea_int_precip;
  if (!strcmp (name, "int_evap"))
    return airsea_int_evap;
  if (!strcmp (name, "int_fwf"))
    return airsea_int_fwf;
  if (!strcmp (name, "int_swr"))
    return airsea_int_swr;
  if (!strcmp (name, "int_heat"))
    return airsea_int_heat;
  if (!strcmp (name, "int_total"))
    return airsea_int_total;
  if (!strcmp (name, "bio_drag_scale"))
    return airsea_bio_drag_scale;
  if (!strcmp (name, "bio_albedo"))
    return airsea_bio_albedo;
  if (!strcmp (name, "swr_method"))
    return airsea_swr_method;
  if (!strcmp (name, "albedo_method"))
    return airsea_albedo_method;
  if (!strcmp (name, "const_albedo"))
    return airsea_const_albedo;
  if (!strcmp (name, "fluxes_method"))
    return airsea_fluxes_method;
  if (!strcmp (name, "back_radiation_method"))
    return airsea_back_radiation_method;
  if (!strcmp (name, "rain_impact"))
    return airsea_rain_impact;
  if (!strcmp (name, "calc_evaporation"))
    return airsea_calc_evaporation;
  return HUGE;
}
