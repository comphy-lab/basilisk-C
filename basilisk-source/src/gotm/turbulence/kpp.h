// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/turbulence/kpp.F90
extern realtype __kpp_MOD_zsbl;
#define kpp_zsbl __kpp_MOD_zsbl
extern realtype __kpp_MOD_zbbl;
#define kpp_zbbl __kpp_MOD_zbbl
extern realtype __kpp_MOD_ric;
#define kpp_ric __kpp_MOD_ric
extern logical __kpp_MOD_kpp_sbl;
#define kpp_kpp_sbl __kpp_MOD_kpp_sbl
extern logical __kpp_MOD_kpp_bbl;
#define kpp_kpp_bbl __kpp_MOD_kpp_bbl
extern logical __kpp_MOD_kpp_interior;
#define kpp_kpp_interior __kpp_MOD_kpp_interior
extern logical __kpp_MOD_clip_mld;
#define kpp_clip_mld __kpp_MOD_clip_mld

extern void __kpp_MOD_init_kpp (
  integer * namlst,
  char * fn,
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * kpp_g,
  realtype * kpp_rho_0,
  long int _slfn
);
static inline void kpp_init_kpp (
  integer * namlst,
  char * fn,
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * kpp_g,
  realtype * kpp_rho_0) {
  __kpp_MOD_init_kpp (namlst, fn, nlev, h0, h, kpp_g, kpp_rho_0, strlencheck(fn));
}

extern void __kpp_MOD_do_kpp (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * nnt,
  realtype * nns,
  realtype * ss,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f
);
static inline void kpp_do_kpp (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * nnt,
  realtype * nns,
  realtype * ss,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f) {
  __kpp_MOD_do_kpp (nlev, h0, h, rho, u, v, nn, nnt, nns, ss, u_taus, u_taub, tflux, btflux, sflux, bsflux, trad, brad, f);
}

extern void __kpp_MOD_interior (
  integer * nlev,
  realtype * nn,
  realtype * nnt,
  realtype * nns,
  realtype * ss
);
static inline void kpp_interior (
  integer * nlev,
  realtype * nn,
  realtype * nnt,
  realtype * nns,
  realtype * ss) {
  __kpp_MOD_interior (nlev, nn, nnt, nns, ss);
}

extern void __kpp_MOD_surface_layer (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f
);
static inline void kpp_surface_layer (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f) {
  __kpp_MOD_surface_layer (nlev, h0, h, rho, u, v, nn, u_taus, u_taub, tflux, btflux, sflux, bsflux, trad, brad, f);
}

extern void __kpp_MOD_bottom_layer (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f
);
static inline void kpp_bottom_layer (
  integer * nlev,
  realtype * h0,
  realtype * h,
  realtype * rho,
  realtype * u,
  realtype * v,
  realtype * nn,
  realtype * u_taus,
  realtype * u_taub,
  realtype * tflux,
  realtype * btflux,
  realtype * sflux,
  realtype * bsflux,
  realtype * trad,
  realtype * brad,
  realtype * f) {
  __kpp_MOD_bottom_layer (nlev, h0, h, rho, u, v, nn, u_taus, u_taub, tflux, btflux, sflux, bsflux, trad, brad, f);
}

extern void __kpp_MOD_wscale (
  realtype * bfsfc,
  realtype * u_taus,
  realtype * d,
  realtype * wm,
  realtype * ws
);
static inline void kpp_wscale (
  realtype * bfsfc,
  realtype * u_taus,
  realtype * d,
  realtype * wm,
  realtype * ws) {
  __kpp_MOD_wscale (bfsfc, u_taus, d, wm, ws);
}

extern void __kpp_MOD_clean_kpp (

);
static inline void kpp_clean_kpp (void) {
  __kpp_MOD_clean_kpp();
}
realtype kpp_get_global (const char * name) {
  if (!strcmp (name, "zsbl"))
    return kpp_zsbl;
  if (!strcmp (name, "zbbl"))
    return kpp_zbbl;
  if (!strcmp (name, "ric"))
    return kpp_ric;
  if (!strcmp (name, "kpp_sbl"))
    return kpp_kpp_sbl;
  if (!strcmp (name, "kpp_bbl"))
    return kpp_kpp_bbl;
  if (!strcmp (name, "kpp_interior"))
    return kpp_kpp_interior;
  if (!strcmp (name, "clip_mld"))
    return kpp_clip_mld;
  return HUGE;
}
