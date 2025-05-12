void print_help(){
  fprintf(stderr, "Usage: ./kelvin_parallel <LEVEL> <GEOMETRYFILE> <VOXELSIZE> <STRUTDIAM> <REYNOLDSNUM> <STOKESNUM>\n");
  fprintf(stderr, "\t<LEVEL>\t\t\tdefines resolution=2^(LEVEL)\n");
  fprintf(stderr, "\t<GEOMETRYFILE>\t\tset binary file for Cacluclation: has to fit with LEVEL\n");  
  fprintf(stderr, "\t<STRUTDIAM>\t\tset physical strut diam (grid/geometry dependend)\n");
  fprintf(stderr, "\t<REYNOLDSNUM>\t\tset set Reynolds number\n");  
  fprintf(stderr, "\n\tTo Read Config from file:  mpirun -np 8 ./kelvin_parallel  $(cat config.cfg | tail -n +LINE | head -n 1)");
  exit(1);
}

/** Use Multigrid and the corresponding exporter*/
#include "grid/multigrid3D.h"
#include "output_vti_pv590.h"

/** Use embed to describe the complex boundary geometry*/
#include "embed.h"
/** Use the incompressible Navier-Stokes solver - defined at cell center*/
#include "navier-stokes/centered.h"
/** enable performance monitoring*/
#include "navier-stokes/perfs.h"
/** enables controlled ending of the simulation befor slurm kills it*/
#include "maxruntime.h"
/** Read binary distance field*/
#include "readxyz.h" // read_double_data()
/** function for embed force caluclation */
#include "../../ghigo/src/myembed-force.h" 

// stores avg velocity of the previous time step
double ux_previous_savg = 0.;
// defines time step when simulation will end, can be overwritten
double t_end = 10000.;

/** variables red from config file */
int LEVEL;
char geometryFile[160];

/**  Macro for time averaging
 *   x = avg_field, y = field to average  */
#define TIME_AVG(x,y,t_res) ((x)=((x)*(t-(t_res))/(t-(t_res)+dt) + (y)*dt/(t-(t_res)+dt)))

/** Reynolds stress tensor of current time step*/
vector ReynoldsStressX[];
vector ReynoldsStressY[];
vector ReynoldsStressZ[];

/** time averaged Reynolds stress tensor */
vector ReynoldsStressX_tavg[];
vector ReynoldsStressY_tavg[];
vector ReynoldsStressZ_tavg[];
/** Distance field (wall distance)*/
scalar d[];
/** turbuletic kinetic energy of current time step*/
scalar k_t[];
/** time averaged turbuletic kinetic energy */
scalar k_t_tavg[];

/** x-direction of velocity of the last time step*/
scalar u_last[];
/** velocity field of the last time step*/
vector u_prev[];
/** velocity field of the last time step*/
vector u_tavg[];

/** Vector for Forces on Surface*/
vector Fp_v[];
vector Fmu_v[];
vector Fp_v_tavg[];
vector Fmu_v_tavg[];

/** Average velocity */
double muv_stavg = 0.;
/** Value Reset time */
double t_reset = 0.;
double t_reset_muv = 0.;
/** Defines the size of the domain in physical units, set via command line arguments */
double STRUTDIAM;
double Re_soll;
//Verweilzeit
double t_comp = HUGE;
double t_start = 0.;
double t_start0 = 0.;
double muv_start = 1.; // some value not to zero
double muv_stavg_total = 0.;
double muv_stavg_interval = 0.;
double A_fluk = 0.;


//~ int once = 0;

/** Global face vector field for the variable */
face vector muv[];

/** Main function */
int main(int argc, char * argv[])
{  
  /**  Check, if maximum runtime is given as an command line arg (-m "HH:MM:SS") */
  maxruntime (&argc, argv);
  /**  Check for command line arg*/
  if (argc >= 5) {
    LEVEL = atoi(argv[1]);
    sprintf(geometryFile,"%s",argv[2]);    
    STRUTDIAM = atof(argv[3]);
    Re_soll = atof(argv[4]);  
    assert(LEVEL > 4);
    assert(Re_soll > 1e-6);    
    fprintf(stderr,\
      "LEVEL: %i GEOMETRYFILE:%s STRUCTDIAM:%g REYNOLDSNUM:%g\n",\
      LEVEL, geometryFile, STRUTDIAM, Re_soll);
  } else {    
    print_help();
    exit(1);
  }
  /** check if npe() is 8^[1..2..3..] to ensure correct domain shape*/
  if (!(powl((int)round(cbrt(npe())),3) == npe())) { //       
    fprintf(stdout,"Richtige Anzahl an Prozessen angeben\n");
    fprintf(stdout,"Usage: mpirun -np X ... , X = 8^[0..1..2..3..4...]\n");    
    exit(1);
  }
    //~ system("mkdir -p /media/data/vti/");
    //~ system("mkdir -p vti/");
    fprintf(stdout,\
      "Prozess:%04i: Erzeuge FractionField aus DistanceField für MG auf Level %i\n"\
      ,pid(), LEVEL);
    
  /** Set Resolution for the longest axis, 2^LEVEL*/
  init_grid (1<<(LEVEL));
  /** Set Domain size*/  
  size (1<<(LEVEL));
  
  /**
  * Periodische Randbedingung in allen Raumrichtungen
  */
  foreach_dimension()
    periodic(right);
 
  /**
   *  Set max. Timestep, (if needed) */
    const face vector g[] = {1.,0.,0.};
    a = g;
    mu = muv;  
  TOLERANCE=1e-5;
  /** Limit max. Timstep */
  DT=1e-1; 
  
  run();
}


/** Set BC on the embed filter surface
 *  n: normal
 *  t: tangential
 *   r: tangential_2 (for 3D)
 * */
    
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
#if dimension==3
u.r[embed] = dirichlet(0.);
#endif

/** Initializing */ 
event init (i = 0) {
  //~ CFL=0.3; // Set CFL number if needed
  FILE * fp ;  
    fprintf(stderr,"here\n");
  char fileToRead[160];
    
  /** Open and read Distance field and compute fraction (cs) from it */
  sprintf(fileToRead, "%s", geometryFile);        
  fp = fopen(fileToRead, "rb");
  if (!fp) { fprintf(stdout,"DistanceField not found\n"); return 1; } 
  
  read_double_data(fp, d, LEVEL);          
  fclose(fp);
  
  /** Calculate Distancefield on vertices (average)*/
  vertex scalar phi[];
  foreach_vertex(){
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
    d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  }
  boundary ({phi});
    
  /** Calculate Fractions of the embed geometry) */
  fractions (phi,cs,fs);       
  fractions_cleanup(cs,fs);
  boundary({cs,fs}); // probably not needed
  restriction({cs,fs}); // probably not needed
     
  foreach(){    
    u.x[] = 1.*cs[];
  }
  /** Output STL-Geometry */
  //stl_output_binary (cs, "csg.stl"); // broken at the moment
}

event properties (i++)
{
  //~ if (once == 0)
  {
    double ux_savg = normf(u.x).avg;
    
    //~ for (int ii=1; ii<8; ++ii){
      //~ u_hist[ii] = u_hist[ii-1];      
    //~ }
    //~ u_hist[0] = ux_savg;
    //~ double sum = 0.;
    //~ for (int ii=0; ii<8; ++ii){
      //~ sum += u_hist[ii];
    //~ }
    //~ double ux_avg = sum/8.;
    
    //~ double grad = -(ux_savg -ux_previous_savg)/dt;
    
    foreach_face(){
       //~ muv.x[] = fm.x[]*STRUTDIAM*(ux_savg+(1-exp(-30*i/10))*grad*dt)/Re_soll;
       //~ muv.x[] = fm.x[]*STRUTDIAM*(ux_savg+ux_previous_savg+(1.-exp(-5.*i/10.))*grad*dt)/(2*Re_soll);
       //~ muv.x[] = fm.x[]*STRUTDIAM*ux_avg/Re_soll;
       //~ muv.x[] = fm.x[]*STRUTDIAM*(ux_savg+ux_previous_savg)/(2*Re_soll);  
       muv.x[] = fm.x[]*STRUTDIAM*ux_savg/Re_soll;       
    }
  ux_previous_savg = ux_savg;        
}
}

//~ event avg_velocity(i++)
//~ {
  //~ foreach()
    //~ foreach_dimension()
      //~ TIME_AVG(u_tavg.x[],u.x[],t_reset);
  //~ boundary((scalar*){u_tavg});  
//~ }

event calc_turbulent_energy(i++)
{
  foreach(){
    foreach_dimension() {
      TIME_AVG(u_tavg.x[],u.x[],t_reset);
      ReynoldsStressX.x[] = sq(u.x[] - u_tavg.x[]);
    }
    ReynoldsStressX.y[] = (u.x[] - u_tavg.x[]) * (u.y[] - u_tavg.y[]);
    ReynoldsStressY.x[] = ReynoldsStressX.y[];
    
    ReynoldsStressX.z[] = (u.x[] - u_tavg.x[]) * (u.z[] - u_tavg.z[]);
    ReynoldsStressZ.x[] = ReynoldsStressX.z[];
    
    ReynoldsStressY.z[] = (u.y[] - u_tavg.y[]) * (u.z[] - u_tavg.z[]);
    ReynoldsStressZ.y[] = ReynoldsStressY.z[];
    
    k_t[] = 1./2.*(sq(ReynoldsStressX.x[]) \
                 + sq(ReynoldsStressY.y[]) \
                 + sq(ReynoldsStressZ.z[]));    
    
    TIME_AVG(k_t_tavg[],k_t[],t_reset);
    /** Calculate the average Reynolds stress tensor*/
    foreach_dimension(){
      TIME_AVG(ReynoldsStressX_tavg.x[],ReynoldsStressX.x[],t_reset);
      TIME_AVG(ReynoldsStressY_tavg.x[],ReynoldsStressY.x[],t_reset);
      TIME_AVG(ReynoldsStressZ_tavg.x[],ReynoldsStressZ.x[],t_reset);  
  
    }
  }
}

event calc_force(i++) {
  coord Fp, Fmu;
  double Fp_x = 0., Fp_y = 0., Fp_z = 0., Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;
  foreach (reduction(+:Fp_x) reduction(+:Fp_y) reduction(+:Fp_z)
     reduction(+:Fmu_x) reduction(+:Fmu_y) reduction(+:Fmu_z)) {
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn = area*embed_interpolate_3D (point, p, b);
      foreach_dimension(){
        double tmp = Fn*n.x;    
        Fp_v.x[] = tmp;
        TIME_AVG(Fp_v_tavg.x[],tmp,t_reset);          
        Fp_x += tmp;
      }
      if (constant(mu.x) != 0.) {
        double mua = 0., fa = 0.;
        foreach_dimension() {
          mua += mu.x[] + mu.x[1];
          fa  += fs.x[] + fs.x[1];
        }
        mua /= fa;

        coord dudn = embed_gradient (point, u, b, n);
  #if dimension == 2
        foreach_dimension()
          double tmp = area*mua*(dudn.x*(sq (n.x) + 1.) +
                       dudn.y*n.x*n.y);
          Fmu_v.x[] = tmp;
          TIME_AVG(Fmu_v_tavg.x[],tmp,t_reset);
          Fmu_x    -= tmp;           
  #else // dimension == 3
        foreach_dimension(){
          double tmp = area*mua*(dudn.x*(sq (n.x) + 1.) +
                       dudn.y*n.x*n.y +
                       dudn.z*n.x*n.z);
          Fmu_v.x[]  = tmp;
          TIME_AVG(Fmu_v_tavg.x[],tmp,t_reset);
          Fmu_x    -= tmp;
        }
  #endif // dimension
      }
    }
  }
  foreach_dimension() {
    Fp.x = Fp_x;
    Fmu.x = Fmu_x;
  }
}
    
/**
 * Laufzeitstatistiken
 */
event logfile (i++){  
  /**
   * Print relevant information on screen
   * - Filtration efficiency: np_filtered/np
   * - Pressure Gradient: acc
   * - normf(u.x).avg
   * - Re
   */ 
  double ux_savg = normf(u.x).avg;
  double muv_savg = normf(muv.x).avg; // avg over muv.x 
  
  double RE = STRUTDIAM*ux_savg/(muv_savg+SEPS);  
  TIME_AVG(muv_stavg,muv_savg,t_reset_muv);
  if (i == 10){
    t_reset_muv = t;
    muv_savg = 0.;
  }
  
  if (i==0) fprintf(stderr,"%-10s %-10s %-10s %-10s %-10s %-10s\n",
                            "Iteration","Time","u.x.avg","muv.x.avg",\
                            "Re","rel. change muv.x");    
  fprintf(stderr,"%10d %-10.6f %-10.6g %-10.6f %-10.6f %-10.6g\n",\
  i,t,ux_savg,muv_savg,RE,fabs((muv_savg-muv_stavg)/muv_stavg));
  
  #if 0
  /** Check if the viskosity is somewhat converged */
  if ((once == 0) && (fabs((muv_savg-muv_stavg)/muv_stavg) < 5e-4) && (i>300)) {
    once = 1;
    fprintf(stdout,"#############################\n");
    fprintf(stdout,"*  Reset    i:%i, t:%f      \n",i,t);
    fprintf(stdout,"#############################\n");
    foreach(){ // Reset time averaged values
      foreach_dimension(){
        ReynoldsStressX_tavg.x[] = 0.;
        ReynoldsStressY_tavg.x[] = 0.;
        ReynoldsStressZ_tavg.x[] = 0.;    
      }      
      k_t_tavg[] = 0.;      
      u_tavg.x[] = 0.;
    }
    /** caluculate 5 Verweilzeiten */
    double t_stay += L0 / (normf(u.x).avg + SEPS);    
    t_end = t + 5.*t_stay;  
    t_reset = t;
  }
  #endif
  /** Check if case convergef */
  if (i==300){
    muv_stavg_total = 0.;
    t_start0 = t;
  }
  if (i==350){
    t_start = t;    
    t_comp = L0 / (normf(u.x).avg) + t_start;
    muv_start = muv_savg;
    muv_stavg_interval = 0.;    
    A_fluk = 0.;
    fprintf(stdout,"Next Comparison Interval: %g\n",t_comp);
  }
  TIME_AVG(muv_stavg_total,muv_savg,t_start0); // intial averaged values
  TIME_AVG(muv_stavg_interval,muv_savg,t_start); // interval averaged values
  A_fluk += fabs(muv_stavg_interval - muv_stavg_total)/muv_stavg_total * dt;
  
  if (t > t_comp){
    fprintf(stdout,"New Comparison Interval\n");    
    double A_ref = (t - t_start)*muv_start;
    fprintf(stdout,"A_fluk: %g ,A_ref: %g\n",A_fluk,A_ref);
    fprintf(stdout,"New Comparison Interval: A_fluk/A_ref: %g\n",A_fluk/A_ref);
    if (A_fluk/A_ref < 1e-3){
      fprintf(stdout,"Simulation konvergiert!\n");
      return 1;
    }
    A_fluk = 0.;    
    
    t_start = t;
    t_comp = L0 / (normf(u.x).avg) + t_start; // increase t_comp by 1 verweilzeit
    muv_stavg_interval = 0.;
    muv_start = muv_savg;
    fprintf(stdout,"Next Comparison Interval: %g\n",t_comp);
    fflush(stdout);
  }
}
/**
 * Dump
 * use 'qcc -D_DUMP ...' to activate
 */
#if _DUMP
event dumping (i+=100)
{
  char name[80];
  sprintf (name, "datadump-%d", i);
  dump(file = name);  
}
#endif
#if 1
//~ event snapshot (i+=100) 
event snapshot (i+=20) 
{    
  scalar mpi_proc[], l_ref[];
  
  foreach_cell(){
    l_ref[] = level;    
    mpi_proc[] = pid();    
  }
  vector U[];
  foreach()
    foreach_dimension()
      U.x[] = cs[]*u.x[];
  boundary(all);
        
  char path[80];
  sprintf(path, "vti_level%i_Re%i_npr%i", LEVEL,(int)Re_soll,npe());  
  if (i == 0)
  {
    char mkdir_cmd[180];
    sprintf(mkdir_cmd, "mkdir -p %s", path);
    system(mkdir_cmd);  
    sleep(2); // wait for folder to be created
  }
  char prefix[80];
  sprintf(prefix, "data_%06d", i);      
  output_vti((scalar *){cs, p, mpi_proc,d,k_t,k_t_tavg}, \
         (vector *){U,u_tavg,Fp_v,Fmu_v,Fp_v_tavg,Fmu_v_tavg,ReynoldsStressX,ReynoldsStressY,ReynoldsStressZ,ReynoldsStressX_tavg,ReynoldsStressY_tavg,ReynoldsStressZ_tavg}, path, prefix, i, t);    
}
#endif

event ending(i++)
  if (t >= t_end) { // Check, wether the Verweilzeiten are over
    {
      scalar mpi_proc[], l_ref[];
  
      foreach_cell(){
        l_ref[] = level;    
        mpi_proc[] = pid();    
      }
      vector U[];
      foreach()
        foreach_dimension()
          U.x[] = cs[]*u.x[];
      boundary(all);
          
      //~ // char path[]="vti/";
      char path[80];
      sprintf(path, "vti_level%i_Re%i_npr%i", LEVEL,(int)Re_soll,npe());  
      if (i == 0)
      {
        char mkdir_cmd[180];
        sprintf(mkdir_cmd, "mkdir -p %s", path);
        system(mkdir_cmd);  
        sleep(3); // wait for folder to be created
      }
      char prefix[80];
      sprintf(prefix, "data_%06d", i);      
      output_vti((scalar *){cs, p, mpi_proc,d,k_t,k_t_tavg}, \
             (vector *){U,u_tavg,Fp_v,Fmu_v,Fp_v_tavg,Fmu_v_tavg,ReynoldsStressX,ReynoldsStressY,ReynoldsStressZ,ReynoldsStressX_tavg,ReynoldsStressY_tavg,ReynoldsStressZ_tavg}, path, prefix, i, t);      
    }    
    printf("Endzeit erreicht, Simulation beendet\n");
    return 1;
  }
  
event simulation_ending(t=10000)
{ // Simulation will end at t=100 (large placeholder value, otherwise simulation will stop after 10 Steps
  printf("Simulation beendet\n");
  return 1;
}

