// Generated automatically by gotm.awk from
// /home/popinet/local/src/GOTM-5.2.1/src/util/time.F90
extern char * __time_MOD_timestr;
#define time_timestr __time_MOD_timestr
extern char * __time_MOD_start;
#define time_start __time_MOD_start
extern char * __time_MOD_stop;
#define time_stop __time_MOD_stop
extern realtype __time_MOD_timestep;
#define time_timestep __time_MOD_timestep
extern realtype __time_MOD_fsecs;
#define time_fsecs __time_MOD_fsecs
extern realtype __time_MOD_simtime;
#define time_simtime __time_MOD_simtime
extern realtype __time_MOD_fsecondsofday;
#define time_fsecondsofday __time_MOD_fsecondsofday
extern integer __time_MOD_julianday;
#define time_julianday __time_MOD_julianday
extern integer __time_MOD_secondsofday;
#define time_secondsofday __time_MOD_secondsofday
extern integer __time_MOD_yearday;
#define time_yearday __time_MOD_yearday
extern integer __time_MOD_timefmt;
#define time_timefmt __time_MOD_timefmt
// static const integer time_timestepkind = selected_int_kind;
extern integer __time_MOD_minn;
#define time_minn __time_MOD_minn
extern integer __time_MOD_maxn;
#define time_maxn __time_MOD_maxn
extern integer __time_MOD_jul2;
#define time_jul2 __time_MOD_jul2
extern integer __time_MOD_secs2;
#define time_secs2 __time_MOD_secs2

extern void __time_MOD_init_time (
  timestepkind * minn,
  timestepkind * maxn
);
static inline void time_init_time (
  timestepkind * minn,
  timestepkind * maxn) {
  __time_MOD_init_time (minn, maxn);
}

extern void __time_MOD_calendar_date (
  integer * julian,
  integer * yyyy,
  integer * mm,
  integer * dd
);
static inline void time_calendar_date (
  integer * julian,
  integer * yyyy,
  integer * mm,
  integer * dd) {
  __time_MOD_calendar_date (julian, yyyy, mm, dd);
}

extern void __time_MOD_julian_day (
  integer * yyyy,
  integer * mm,
  integer * dd,
  integer * julian
);
static inline void time_julian_day (
  integer * yyyy,
  integer * mm,
  integer * dd,
  integer * julian) {
  __time_MOD_julian_day (yyyy, mm, dd, julian);
}

extern void __time_MOD_update_time (
  timestepkind * n
);
static inline void time_update_time (
  timestepkind * n) {
  __time_MOD_update_time (n);
}

extern void __time_MOD_read_time_string (
  char * timestr,
  integer * jul,
  integer * secs,
  long int _sltimestr
);
static inline void time_read_time_string (
  char * timestr,
  integer * jul,
  integer * secs) {
  __time_MOD_read_time_string (timestr, jul, secs, strlencheck(timestr));
}

extern void __time_MOD_write_time_string (
  integer * jul,
  integer * secs,
  char * timestr,
  long int _sltimestr
);
static inline void time_write_time_string (
  integer * jul,
  integer * secs,
  char * timestr) {
  __time_MOD_write_time_string (jul, secs, timestr, strlencheck(timestr));
}

extern realtype __time_MOD_time_diff (
  integer * jul1,
  integer * secs1,
  integer * jul2,
  integer * secs2
);
static inline realtype time_time_diff (
  integer * jul1,
  integer * secs1,
  integer * jul2,
  integer * secs2) {
  return __time_MOD_time_diff (jul1, secs1, jul2, secs2);
}

extern void __time_MOD_sunrise_sunset (
  realtype * latitude,
  realtype * declination,
  realtype * sunrise,
  realtype * sunset
);
static inline void time_sunrise_sunset (
  realtype * latitude,
  realtype * declination,
  realtype * sunrise,
  realtype * sunset) {
  __time_MOD_sunrise_sunset (latitude, declination, sunrise, sunset);
}

extern void __time_MOD_print_state_time (

);
static inline void time_print_state_time (void) {
  __time_MOD_print_state_time();
}
realtype time_get_global (const char * name) {
  if (!strcmp (name, "timestep"))
    return time_timestep;
  if (!strcmp (name, "fsecs"))
    return time_fsecs;
  if (!strcmp (name, "simtime"))
    return time_simtime;
  if (!strcmp (name, "fsecondsofday"))
    return time_fsecondsofday;
  if (!strcmp (name, "julianday"))
    return time_julianday;
  if (!strcmp (name, "secondsofday"))
    return time_secondsofday;
  if (!strcmp (name, "yearday"))
    return time_yearday;
  if (!strcmp (name, "timefmt"))
    return time_timefmt;
  if (!strcmp (name, "minn"))
    return time_minn;
  if (!strcmp (name, "maxn"))
    return time_maxn;
  if (!strcmp (name, "jul2"))
    return time_jul2;
  if (!strcmp (name, "secs2"))
    return time_secs2;
  return HUGE;
}
