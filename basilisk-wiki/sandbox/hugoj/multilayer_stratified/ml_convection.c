/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

  * Compilation and run (with mpirun)
  make

  * Compilation and run (gpu)
  __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_convection.gpu.tst

  * HPC source file generation
  make _ml_convection.c

  * Running on HPC

  -> first compile on hpc with
  mpicc -std=c99 -O2 _*.c

  -> then run (in a slurm file)
  srun ./ml_convection

*/


#include "grid/multigrid.h"
//#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/diffusion.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

//#define g_ 9.81

/*

TO DO : 
- add dimensions so that they can be saved in the netcdf


- aspect ratio: H/L with H height of OBL/horizontal size of plume
- Convective Richardson: N**2.H**2/(w*)**2
- convective scale w*=(-B0.h)**(1/3)
- surface buoyancy flux B0=g*betaT*Q/(rho.cp)

- convective time t* = H/w*
- analytical law: h(t) = sqrt(2.8/N**2 . -B0.t) (Van Roeckel)

*/

/*
DEFAULT PARAMETERS

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";          // file name of output
// -> Initial conditions
double strat = 0.000002 [0,-2];       // [s-2] N^2 stratification
double Ts = 20. [0,0,1];              // [K] Surface temperature (arbitrary)
// -> Forcing
double qt = 100. [-2,-1,0,1];         // [W.m-2] Heat flux
// -> Domain definition
int N_grid = 5 [0];                   // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];                 // domain size
int N_layer = 2 [0];                  // number of layers
double h0 = 1.0 [1];                  // depth of water
// -> Runtime parameters
int restart = 0 [0];                  // 1: restart, 0: no restart
double tend = 2.0 [0,1];              // end time of simulation
// -> saving outputs
double dtout = 2.0 [0,1];             // dt for output in netcdf
int pad = 4 [0];                      // number of 0-padding for ouput files
int nout = 1 [0];                     // number of the outfile
char fileout[100];                    // name of outfile
// -> physical properties
double nu0 = 0.00025 [1,-2,0,0,1];    // Viscosity for vertical diffusion
double thetaH = 0.5 [0];              // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025. [-3,0,0,0,1];     // [kg.m-3] reference density
double cp = 4.2e3 [0,0,-1,1,-1];      // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4 [0,0,-1];         // [K-1] Thermal expansion coeff for water
double D = 1.5e-5 [2,-1];             // [m2.s-1] Scalar vertical diffusion coeff
double T0 = 20. [0,0,1];              // [°C] Reference temperature
double Trand = 0.1 [0,0,1];           // [°C] Random temperature perturbution
const double g_ = 9.81 [1,-2];        // [m.s-2] Gravity
const double nodim = 1. [0];
#define drho(T) (betaT*(T0-T))        // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
#define Tini(z) strat/(g_*betaT)*z + Ts

// double drho0 = 1e-3, T0 = 13.25, T1 = 8., H0 = 100.;
// #define drho(T) (drho0*(T0 - T)/(T0 - T1))
// #define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
//#define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
#include "layered/dr.h"

// diag
double *T_profile;
double *W_profile;
double *Tflx_profile;
double dt_mean = 1.;

static FILE * fp1;
static FILE * fp2;



int main(int argc, char *argv[])  
{
 
  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
  add_param("dtout", &dtout, "double");
  add_param("strat", &strat, "double");
  add_param("Ts", &Ts, "double");
  add_param("rho0", &rho0, "double");
  add_param("cp", &cp, "double");
  add_param("betaT", &betaT, "double");
  add_param("qt", &qt, "double");
  add_param("dt_mean", &dt_mean, "double");
  
  show_dimension (u);
  dimensional (T_profile, T0);

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  // Settings solver values from namlist values
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 
  //TOLERANCE = 1e-8;
  //CFL_H = 40; 
  //CFL=0.5;

  // Boundary condition
  origin (-L0/2., -L0/2.);
  
  //periodic (right);
  
  periodic (top);
  periodic (left);

  // T profile
  T_profile = (double*)calloc(nl, sizeof(double));
  W_profile = (double*)calloc(nl, sizeof(double));
  Tflx_profile = (double*)calloc(nl, sizeof(double));
  fp1  = fopen("T_profile.dat","w"); // reset file
  fclose(fp1);
  fp2  = fopen("Tflx_profile.dat","w"); // reset file
  fclose(fp2);

  fprintf (stderr, "Read in parameters!\n");
  run();
  
}

event init(i =  0) {
  geometric_beta (0., true); // if !=0, varying layer thickness
  // step 1: set eta and h
  foreach() {
    zb[] = -h0;
    eta[] = 0.;
    double H = - zb[];
    
    foreach_layer() {
      h[] = H/nl;
    } 
  }
  
  // step 2: remap
  vertical_remapping (h, tracers);
  // step 3: set currents
  foreach() {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      u.x[] = 0.;
      u.y[] = 0.;
      w[] = 0.;
      T[] = Tini(z) - Trand + (rand() / (RAND_MAX+1.0) * 2 * Trand) ;
      //T[] = T0(z);
      z += h[]/2.;
    }
  }
  #if dimension==1
    create_nc((scalar *){zb, h, u.x, w, eta, T}, file_out);
  #else
    create_nc((scalar *){zb, h, u.x, u.y, w, eta, T}, file_out);
  #endif
  
  fprintf (stderr,"Done initialization!\n");

}

// vertical diffusion on T
// This implements a heat flux at the surface
event viscous_term (i++)
{
  foreach()
    vertical_diffusion (point, h, T, dt, D, qt/(D*rho0*cp), 0., 0.);
    //TODO: change s_b for the value at bottom of stratification, instead of
    //      this 0.°C that is cooling down the fluid.
}

// Writing a 4D netcdf file
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  write_nc();
}



double* h_avg(scalar var, double* profile){
  /*
  This function computes the horizontal average of var.
  
  INPUTS:
    var: scalar (C Basilisk), the variable to average.
    profile: array of double of length nl
  OUTPUTS:
    an array of double with the average values inside
   */
  foreach(reduction(+:profile[:nl]))
    foreach_layer(){
    #if dimension==1
      profile[point.l] += var[] / N;// (N*N); // * dt / dt_mean;
    #else
      profile[point.l] += var[] / (N*N);
    #endif
    }
  return profile;
}

int write_profile(char* name, double* profile, FILE* fp){
  /*
  docstring

  */
  // main worker is writing the file
  if (pid()==0) {
    fp  = fopen(name,"a");
    if (fp == NULL){
      fprintf(stderr, "Error opening file %s", name);
      return 2;
    }
    for (int i=0; i<nl; ++i) {
      fprintf (fp, "%f %d %g\n", t, i, profile[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
  }
  return 0;
}


// This event compute layer average of T, w
event compute_horizontal_avg (t+=dt_mean; t<=tend+1e-10){
  T_profile = h_avg(T, T_profile);
  W_profile = h_avg(w, W_profile);
}

event compute_flx_wT(t+=dt_mean, t<=tend+1e-10){
  foreach(reduction(+:Tflx_profile[:nl]))
    foreach_layer(){
    #if dimension==1
      Tflx_profile[point.l] += ((w[]-W_profile[point.l])*(T[]-T_profile[point.l])) / N;
    #else
      Tflx_profile[point.l] += ((w[]-W_profile[point.l])*(T[]-T_profile[point.l])) / (N*N);
    #endif
  }
}


// This even writes to a file the layer average
event write_diag(t=0., t+=dt_mean){
    // main worker is writing the file
    // if (pid()==0) {
    //   fp  = fopen("T_profile.dat","a");
    //   if (fp == NULL){
    //     fprintf(stderr, "Error opening file T_profile.dat");
    //     return 2;
    //   }
    //   for (int i=0; i<nl; ++i) {
    //     fprintf (fp, "%f %d %g\n", t, i, T_profile[i]);
    //   }
    //   fprintf(fp,"\n");
    //   fclose(fp);
    // }
    write_profile("T_profile.dat", T_profile, fp1);
    write_profile("Tflx_profile.dat", Tflx_profile, fp2);

    // Reset the profile for all workers
    for (int i=0; i<nl; ++i) {
      T_profile[i] = 0.0;
      Tflx_profile[i] = 0.0;
    }
}


event cleanup(t=end){
  free(T_profile);
  free(Tflx_profile);
}

/**
Results: plots
~~~gnuplot 
# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "T (°C)"
set yrange [0:14]
set xrange [0:200]
set terminal pngcairo size 800,600 enhanced font 'Verdana_r,12'
set output 'T_profile.png'
set size 0.9, 0.9
splot "T_profile.dat" using 1:2:3 with pm3d
unset output
~~~

**/
