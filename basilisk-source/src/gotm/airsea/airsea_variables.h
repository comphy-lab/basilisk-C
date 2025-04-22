// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/airsea_variables.F90
static const realtype airsea_variables_rgas = 287.1;
static const realtype airsea_variables_g = 9.81;
static const realtype airsea_variables_rho_0 = 1025.;
static const realtype airsea_variables_kappa = 0.41;
extern realtype __airsea_variables_MOD_es;
#define airsea_variables_es __airsea_variables_MOD_es
extern realtype __airsea_variables_MOD_ea;
#define airsea_variables_ea __airsea_variables_MOD_ea
extern realtype __airsea_variables_MOD_qs;
#define airsea_variables_qs __airsea_variables_MOD_qs
extern realtype __airsea_variables_MOD_qa;
#define airsea_variables_qa __airsea_variables_MOD_qa
extern realtype __airsea_variables_MOD_l;
#define airsea_variables_l __airsea_variables_MOD_l
extern realtype __airsea_variables_MOD_rhoa;
#define airsea_variables_rhoa __airsea_variables_MOD_rhoa
extern realtype __airsea_variables_MOD_ta;
#define airsea_variables_ta __airsea_variables_MOD_ta
realtype airsea_variables_get_global (const char * name) {
  if (!strcmp (name, "es"))
    return airsea_variables_es;
  if (!strcmp (name, "ea"))
    return airsea_variables_ea;
  if (!strcmp (name, "qs"))
    return airsea_variables_qs;
  if (!strcmp (name, "qa"))
    return airsea_variables_qa;
  if (!strcmp (name, "l"))
    return airsea_variables_l;
  if (!strcmp (name, "rhoa"))
    return airsea_variables_rhoa;
  if (!strcmp (name, "ta"))
    return airsea_variables_ta;
  return HUGE;
}
