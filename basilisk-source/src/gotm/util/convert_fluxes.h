// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/util/convert_fluxes.F90

extern void convert_fluxes_ (
  integer * nlev,
  realtype * g,
  realtype * cp,
  realtype * rho_0,
  realtype * heat,
  realtype * p_e,
  realtype * rad,
  realtype * t,
  realtype * s,
  realtype * tflux,
  realtype * sflux,
  realtype * btflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad
);
static inline void util_convert_fluxes (
  integer * nlev,
  realtype * g,
  realtype * cp,
  realtype * rho_0,
  realtype * heat,
  realtype * p_e,
  realtype * rad,
  realtype * t,
  realtype * s,
  realtype * tflux,
  realtype * sflux,
  realtype * btflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad) {
  convert_fluxes_ (nlev, g, cp, rho_0, heat, p_e, rad, t, s, tflux, sflux, btflux, bsflux, trad, brad);
}
