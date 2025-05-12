/**
# Shoaling solitary wave in the presence of wind

*/
#include "embed.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"


#include "navier-stokes/perfs.h"


#define LOG_FILENAME "log"

#define VIEW_FILENAME_SHOAL "shoal_view.ppm"
#define LEVEL_FILENAME_SHOAL "shoal_levels.ppm"
#define VORT_FILENAME_SHOAL "shoal_vort.ppm"
#define VOF_FILENAME_SHOAL "shoal_vof.ppm"

#define VIEW_FILENAME_WIND "wind_view.ppm"
#define LEVEL_FILENAME_WIND "wind_levels.ppm"
#define VORT_FILENAME_WIND "wind_vort.ppm"
#define VOF_FILENAME_WIND "wind_vof.ppm"

/*
  hotfix 2022-10-18:
  changed init_grid to LEVEL-2 levels from 8
  added minimum error tolerance (since full refine occurs when u*=0)
 */





/** RE = 40000
*/
#define RE 40000


/** LEVEL_SHOAL = 14
*/
#define LEVEL_SHOAL 14


/** LEVEL_WIND = 11
*/
#define LEVEL_WIND 11


/** BO = 4000
*/
#define BO 4000


/** RHOratio = RHOa / RHOw
where RHOw = 1 (1)
      RHOa = 0.000980392156862745 (1/1020)
*/
#define RHOratio 0.000980392156862745


/** MUratio = 1.81e-5/1.00e-3
*/
#define MUratio 0.018099999999999998


/** Y_TOP = # Ha
where Ha = 10 (10)
*/
#define Y_TOP 10


/** A0 = 0.6
*/
#define A0 0.6


/** REstar = 1800
*/
#define REstar 1800


/** Ustar = #REstar * NUa / Ha
where REstar = 1800 (1800)
      Ha = 10 (10)
      NUa = <function  at 0x7f66c2d6dea0> (#MUa / RHOa)
*/
#define Ustar 0.10508754629165151


/** SIM_SIDELEN = 60
*/
#define SIM_SIDELEN 60


/** X0 = 15
*/
#define X0 15


/** TMAX = 1030
*/
#define TMAX 1030


/** FORCING_TERM = # (Ustar**2)/Y_TOP
where Y_TOP = 10 (# Ha)
      Ustar = 0.10508754629165151 (#REstar * NUa / Ha)
*/
#define FORCING_TERM 0.0011043392385599999


/** VORT_MAG_MAX = 10
*/
#define VORT_MAG_MAX 10


/** VIEW_DT_WIND = 1.0
*/
#define VIEW_DT_WIND 1.0


/** DUMP_DT_WIND = 1.0
*/
#define DUMP_DT_WIND 1.0


/** VIEW_DT_SHOAL = 0.05
*/
#define VIEW_DT_SHOAL 0.05


/** DUMP_DT_SHOAL = 0.05
*/
#define DUMP_DT_SHOAL 0.05


/** SOLIT_C = #(G*(H0 + A0))**0.5
where G = 1 (1)
      H0 = 1 (1)
      A0 = 0.6 (0.6)
*/
#define SOLIT_C 1.2649110640673518


/** T_TRANSITION = 1001
*/
#define T_TRANSITION 1001


/** SHOAL_FOCUS_START = T_TRANSITION + 18
where T_TRANSITION = 1001 (1001)
*/
#define SHOAL_FOCUS_START 1019


/** SHOAL_FOCUS_DT = 0.01
*/
#define SHOAL_FOCUS_DT 0.01


/** SHOAL_FOCUS_END = T_TRANSITION + 23
where T_TRANSITION = 1001 (1001)
*/
#define SHOAL_FOCUS_END 1024


/** NUa = #MUa / RHOa
where RHOa = 0.000980392156862745 (1/1020)
      MUa = 5.723722564904766e-07 (#MUw * MUratio)
*/
#define NUa 0.0005838197016202861


/** FLOOR_BOUNDARY_EMBED = f"y+max({SANDBAR_DEPTH},1-max(0,{BEACH_SLOPE}*(x-{FLOOR_WIDTH})))"
where FLOOR_WIDTH = 30 (30)
      SANDBAR_DEPTH = 0.37109375 (0.95/2.56)
      BEACH_SLOPE = 0.0693 (0.0693)
*/
#define FLOOR_BOUNDARY_EMBED y+max(0.37109375,1-max(0,0.0693*(x-30)))


/** DO_ADAPTIVE = ""
*/
#define DO_ADAPTIVE 


/** DO_BODY_FORCE = ""
*/
#define DO_BODY_FORCE 


/** REVERSE_WIND_DIRECTION = ""
*/
#define REVERSE_WIND_DIRECTION 


/** PRESSURE_LEFT_DIRICHLET = ""
*/
#define PRESSURE_LEFT_DIRICHLET 


/** NO_SHOAL_BODYFORCE = ""
*/
#define NO_SHOAL_BODYFORCE 

#define LEVEL ((t < T_TRANSITION)? LEVEL_WIND: LEVEL_SHOAL)

#define WATER_STAGE_THRESHOLD 1


#ifdef DO_PASSIVE_TRACER
#include "tracer.h"
scalar dye[];
scalar * tracers = {dye};
#endif

#ifdef DO_BODY_FORCE
face vector airforce[];
#endif


char _timestr[80];
char * timestamp(){
  struct timeval now;
  gettimeofday(&now,NULL);
  struct tm * ltime = localtime(&now.tv_sec);
  strftime(_timestr, 79, "%F %T", ltime);
  return _timestr;
}

int main()
{

  rho1 = 1.;
  rho2 = RHOratio;
  G.y = -1;
  //DT = 1e-3;
  CFL = 0.5;
  mu1 = 1.0/RE;
  mu2 = 1.0/RE * MUratio;
  f.sigma = 1 / BO; //rho1 * g * h0^2 = 1

  size(SIM_SIDELEN);
  //lower left corner of box:
  origin(0.0,Y_TOP-SIM_SIDELEN); //set upper left corner to (0,Y_TOP)
  
#ifdef DO_BODY_FORCE
  a = airforce;
#endif

  //N = 1 << LEVEL;
  //init_grid(1 << LEVEL);
  init_grid(1 << LEVEL-2);
  run();
}

//=========[initial guess on wind profile]=========
#define smoothstep(x) ((1.0+tanh(0.4 * x))/2.0)
#define boundedlog(x) (x > 1 ? log(x): 0)

double wall_law(double y){
  const double inv_kappa = 1.0/0.41, C = 5.0, buffer_dist = 10.0;
  double yp = y * Ustar * rho2 / mu2;
  //return Ustar * ((yp > buffer_dist)? (inv_kappa * log(yp) + C): yp);
#ifdef REVERSE_WIND_DIRECTION
  return -Ustar * (yp + smoothstep(yp - buffer_dist)*
      (inv_kappa * boundedlog(yp) + C - yp));
#else
  return Ustar * (yp + smoothstep(yp - buffer_dist)*
      (inv_kappa * boundedlog(yp) + C - yp));
#endif
}

//==GN soliton equation, code from MOSTERT & DEIKE (2020):
// http://basilisk.fr/sandbox/wmostert/shallow.c
double waveGN( double x, double y, double a0 ){
  // Try "analytical solution of G-N equations":
  double k = sqrt(3.*a0)/(2.*sqrt(1. + a0));
  return a0*sq(1.0/cosh(k*x)) - y;
}
double detax( double x, double a0 ){
  double tmp1 = -sqrt(3.0)*a0*sqrt(a0)/sqrt(1.0 + a0);
  double tmp2 = sqrt(3.0*a0)/(2.0*sqrt(1.0 + a0));
  double tmp3 = sq(1.0/cosh(tmp2*x));
  double tmp4 = tanh(tmp2*x);
  double deta = tmp1 * tmp3 * tmp4;
  return deta;
}

//=========[BCs: steady inflow/outflow, dirichlet pressure on outflow]=========

#ifdef PRESSURE_LEFT_DIRICHLET
p[right] = neumann(0.0);
pf[right] = neumann(0.0);

p[left] = dirichlet(0.0);
pf[left] = dirichlet(0.0);
#else
p[left] = neumann(0.0);
pf[left] = neumann(0.0);

p[right] = dirichlet(0.0);
pf[right] = dirichlet(0.0);
#endif
u.n[right] = neumann(0.0);
u.n[left] = neumann(0.0);

u.n[top] = dirichlet(0.0);

#ifdef DO_TOP_VEL_CONSTRAINT
u.t[top] = dirichlet(wall_law(Y_TOP));
#else
u.t[top] = neumann(0.0);
#endif

u.n[bottom] = dirichlet(0);

#define SOLIT_ETA(x) waveGN(x-X0,0,A0)
double solit_x_vel(double x){
  double eta = SOLIT_ETA(x);
  return eta / (1.0 + eta);
}
double solit_y_vel(double x, double y){
  double eta = SOLIT_ETA(x);
  double deta = detax(x-X0,A0);
  return -(y + 1.0)*deta/(eta + 1.0) * (1.0 - eta/(1.0 + eta));
}

//bathymetry: no slip (this will be overwritten based on simulation)
int stage;

u.x[embed] = dirichlet(stage == 0? SOLIT_C*(solit_x_vel(x) - 1.0) : 0);
u.y[embed] = dirichlet(stage == 0? SOLIT_C*solit_y_vel(x,y) : 0);


//put tracer in checkerboard pattern.
//left to right flow has values in {-1, -0.5}
//right to left flow has values in {0.5, 1}
#ifdef DO_PASSIVE_TRACER
#define DYE_HORIZK 1./2
#define DYE_FREQ (DYE_HORIZK * 0.3)
#define DYE_VERTK 1./2
dye[left] = (fmod(t * DYE_FREQ, 1) < 0.5 ? -1: 1)*
    (fmod(y * DYE_VERTK, 1) < 0.5 ? -1: 1)*0.25 - 0.75;
dye[right] = (fmod(t * DYE_FREQ, 1) < 0.5 ? -1: 1)*
    (fmod(y * DYE_VERTK, 1) < 0.5 ? -1: 1)*0.25 + 0.75;
#endif

double water_volume(){
    double water = 0;
    foreach(reduction(+:water)){
        water += f[] * dv();
    }
    return water;
}

event init (i = 0)
{
  if(restore("restart")){
    //bugs insue if we do not re-initialize the solid boundary
    
    //which stage are we on?
    if(water_volume() < WATER_STAGE_THRESHOLD){
      //air phase
      solid(cs,fs, y - SOLIT_ETA(x));
      //u.x[embed] = dirichlet(SOLIT_C*(solit_x_vel(x) - 1.0));
      //u.y[embed] = dirichlet(SOLIT_C*solit_y_vel(x,y));
      stage = 0;
    }else{
      //water phase
      solid(cs,fs, FLOOR_BOUNDARY_EMBED);
      //u.x[embed] = dirichlet(0);
      //u.y[embed] = dirichlet(0);
      stage = 1;
    }
    //we should be good to go
    return 0;
  }

  //no restart file: create a new simulation;
  //for now, we only initialize wind simulations

  printf("Failed to load wind_IC. Initializing wind.");
  //take as in wind simulation initialization
  double beta = sqrt(3*A0/(4*(1+A0)));
#define SOLIT_EXPR y - A0*sq(1.0/cosh(beta * (x - X0)))
  const double ymax = A0*1.05;
  // vertical refinement: only above y = 0
  // refine around the wave
  refine(y >= 0 && (
    //refine around wave itself (up to max LEVEL)
    fabs(waveGN(x-X0,y,A0)) < pow(0.5,level-1)
    //refine around the air, reducing one level each distance H0
    //away from the center of the wave (infinity-norm from (X0,0))
    //after two H0 distances
    || level < LEVEL + 2 - floor(max(fabs(x - X0),y))
  )
  && level < LEVEL);
  //embed boundary as soliton
  solid(cs,fs, y - SOLIT_ETA(x));
  //only use air for this stage
  fraction(f,-1);

  foreach(){
#ifdef INIT_ZERO_WIND_AT_WAVE
    u.x[] = (wall_law(y) * smoothstep((y - A0*2)*3.0/0.4)) - SOLIT_C;
#else
    u.x[] = wall_law(y) - SOLIT_C;
#endif
  }
  
}

event switch_stage(t = T_TRANSITION){

  //switch to water stage

  double beta = sqrt(3*A0/(4*(1+A0)));
#define SOLIT_EXPR y - A0*sq(1.0/cosh(beta * (x - X0)))


  refine((
  // refine around the soliton to capture velocity
    (level < LEVEL - floor(fabs(x - X0)*beta) && y > -1.05 && y <= A0)
  // refine the solit boundary
    || (fabs(waveGN(x-X0, y, A0)) < SIM_SIDELEN * pow(0.5,level-1))
  // refine around the bathymetry
    || (fabs(FLOOR_BOUNDARY_EMBED) < SIM_SIDELEN * pow(0.5,level-1)
  )) && level < LEVEL);

  // beach shoaling
  solid(cs,fs, FLOOR_BOUNDARY_EMBED);
  // soliton
  fraction(f,waveGN(x-X0, y, A0));
  //u.x[embed] = dirichlet(0);
  //u.y[embed] = dirichlet(0);
  stage = 1;
  

  //initialize wind via law of the wall, incl. galilean transform
  foreach(){
//==GN soliton equation==
//code from MOSTERT&DEIKE(2020): http://basilisk.fr/sandbox/wmostert/shallow.c
    double eta = waveGN(x-X0, 0, A0);
    double deta = detax(x-X0, A0);
    //u = f * u_solit + (1-f) * u_air
    //u_air = (prior u) + c(+x)
    u.x[] = SOLIT_C * eta / (1.0 + eta) * f[]
	    + (1 - f[]) * (u.x[] + SOLIT_C);
    u.y[] = -(y + 1.0)*SOLIT_C*deta/(eta + 1.0)*(1.0 - eta/(1.0 + eta)) * f[]
	    + (1 - f[]) * u.y[];
  }

}

#ifdef DO_BODY_FORCE
event acceleration(i++){
#ifdef REVERSE_WIND_DIRECTION
  double mag = -FORCING_TERM;//sq(ustar)/(L0-h0)   (air column height)
#else
  double mag = FORCING_TERM;
#endif
#ifdef NO_SHOAL_BODYFORCE
  if(t >= T_TRANSITION){
    mag = 0;
  }
#endif
  foreach_face(x){
#ifdef FORCE_AIR_ONLY
    airforce.x[] += mag*(1.0 - f[]);
#else
    airforce.x[] += mag*(rho2 / rho(f[]));
#endif
  }
}
#endif


#define DIM_VEL_ERR_MAX (1e-2 * max(Ustar,0.04))
#define V_FRAC_ERR_MAX (0.01)

event adapt (i++){
  // Mostert & Deike (2020)
  // ... discretization error of the velocity and VOF fields ...
#ifdef DO_ADAPTIVE
  adapt_wavelet({cs,f,u}, (double[]){0.03,V_FRAC_ERR_MAX,DIM_VEL_ERR_MAX,
      DIM_VEL_ERR_MAX,DIM_VEL_ERR_MAX}, maxlevel = LEVEL);
#endif
  
}
event log_append(i++){
  //KE,PE calculation (Mostert & Dieke), code based on breaking.c
  //code for surface formulation
  double ke = 0.0, pe = 0.0, total_mass = 0.0;
  foreach (reduction(+:ke) reduction(+:pe) reduction(+:total_mass)){
    double vnorm = 0.0;
    foreach_dimension()
      vnorm += sq(u.x[]);
    ke += vnorm * rho(f[]) * dv(); //VOF -> density from two-phase.h
    //assume gravity is only ever in the -y direction
    pe += -rho(f[]) * G.y * y * dv(); 
    total_mass += rho(f[]) * dv();
  }
  ke /= 2;
  
  static FILE * logfp = fopen(LOG_FILENAME, "w");

  fprintf(logfp,"[%s] t = %.5f, ke = %.10e, pe = %.10e, m = %.10e\n",
      timestamp(),t,ke,pe,total_mass);
  fflush(logfp);
}


#define DISP_WATER_VAL (-10)
#define VIEW_BOX {{0,-1.05},{SIM_SIDELEN,Y_TOP}}
event viewer_wind (t = 0; t <= T_TRANSITION; t += VIEW_DT_WIND){
  scalar l[];
  scalar m[];
  scalar vort[];
#ifdef DO_PASSIVE_TRACER
  scalar disp[];
#endif
  foreach(){
    l[] = level;
    m[] = y-waveGN(x-X0,0,A0);
#ifdef DO_PASSIVE_TRACER
    disp[] = (f[] * (DISP_WATER_VAL)) + (1 - f[]) * dye[];
#endif
  }
  vorticity(u,vort);
  //width in pixels (desired height * (width/height))
  // width/height = 1 / (MAXY - MINY)
  int view_width = (int) (256.0 * SIM_SIDELEN / (Y_TOP));
  //print view

#ifdef DO_PASSIVE_TRACER
  static FILE * viewfp = fopen(VIEW_FILENAME_WIND, "w");
  output_ppm(disp,viewfp,min=DISP_WATER_VAL,max=1, linear = true, box = VIEW_BOX,
      n=view_width,mask=m);
#endif
  
  static FILE * voffp = fopen(VOF_FILENAME_WIND, "w");
  output_ppm(f,voffp,min=0,max=1, linear = true, box = VIEW_BOX,
      n=view_width,mask=m);

  static FILE * levelfp = fopen(LEVEL_FILENAME_WIND, "w");

  output_ppm(l,levelfp,min=0,max=LEVEL,n=view_width, box = VIEW_BOX,mask=m);


  static FILE * vortfp = fopen(VORT_FILENAME_WIND, "w");

  output_ppm(vort,vortfp,min=-VORT_MAG_MAX,max=VORT_MAG_MAX,linear=true,
      mask=m,n=view_width,box=VIEW_BOX);
}

event viewer_shoal (t = T_TRANSITION; t += VIEW_DT_SHOAL){
  scalar l[];
  scalar m[];
  scalar vort[];
#ifdef DO_PASSIVE_TRACER
  scalar disp[];
#endif
  foreach(){
    l[] = level;
    m[] = FLOOR_BOUNDARY_EMBED;
#ifdef DO_PASSIVE_TRACER
    disp[] = (f[] * (DISP_WATER_VAL)) + (1 - f[]) * dye[];
#endif
  }
  vorticity(u,vort);
  //width in pixels (desired height * (width/height))
  // width/height = 1 / (MAXY - MINY)
  int view_width = (int) (256.0 * SIM_SIDELEN / (Y_TOP));
  //print view

#ifdef DO_PASSIVE_TRACER
  static FILE * viewfp = fopen(VIEW_FILENAME_SHOAL, "w");
  output_ppm(disp,viewfp,min=DISP_WATER_VAL,max=1, linear = true, box = VIEW_BOX,
      n=view_width,mask=m);
#endif
  
  static FILE * voffp = fopen(VOF_FILENAME_SHOAL, "w");
  output_ppm(f,voffp,min=0,max=1, linear = true, box = VIEW_BOX,
      n=view_width,mask=m);

  static FILE * levelfp = fopen(LEVEL_FILENAME_SHOAL, "w");

  output_ppm(l,levelfp,min=0,max=LEVEL,n=view_width, box = VIEW_BOX,mask=m);


  static FILE * vortfp = fopen(VORT_FILENAME_SHOAL, "w");

  output_ppm(vort,vortfp,min=-VORT_MAG_MAX,max=VORT_MAG_MAX,linear=true,
      mask=m,n=view_width,box=VIEW_BOX);
}

FILE * dump_times_wind = NULL;
int num_dumps_wind = 0;
void do_wind_dump(){
  char fname[80];
  //do dump
  sprintf(fname,"wind_dump/d%04d",num_dumps_wind);
  dump(file = fname);
  
  //store time as (file_index, t) pair.
  if(dump_times_wind == NULL)
    dump_times_wind = fopen("wind_dump/times", "w");
  fprintf(dump_times_wind,"%d %.6f\n",num_dumps_wind,t);

  num_dumps_wind ++;
}
event savestate_wind(t = 0; t < T_TRANSITION; t += DUMP_DT_WIND){
  do_wind_dump();
}

FILE * dump_times_shoal = NULL;
int num_dumps_shoal = 0;
double last_dump_time_shoal = 0;
void do_shoal_dump(){
  if(t - last_dump_time_shoal < 1e-6){
    //don't dump too close to last time
    return;
  }
  char fname[80];
  //do dump
  sprintf(fname,"shoal_dump/d%04d",num_dumps_shoal);
  dump(file = fname);
  
  //store time as (file_index, t) pair.
  if(dump_times_shoal == NULL)
    dump_times_shoal = fopen("shoal_dump/times", "w");
  fprintf(dump_times_shoal,"%d %.6f\n",num_dumps_shoal,t);

  last_dump_time_shoal = t;
  num_dumps_shoal ++;
}
event savestate_shoal(t = T_TRANSITION; t += DUMP_DT_SHOAL){
  do_shoal_dump();
}

#ifdef SHOAL_FOCUS_START
event savestate_shoal_focus(t=SHOAL_FOCUS_START;
        t <= SHOAL_FOCUS_END; t += SHOAL_FOCUS_DT){
  do_shoal_dump();
}
#endif

#if TRACE > 1
event profiling (t += 0.5){
  static FILE * fp = fopen("profiling", "w");
  trace_print(fp,1);
}
#endif

event end (t = TMAX){
  return 1;
}
