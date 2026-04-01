/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

  * Compilation and run (with mpirun)
  make

  * Compilation and run (gpu)
  __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_breaking.gpu.tst

  * HPC source file generation
  make hpc

  * Running on HPC

  -> first compile on hpc with
  mpicc -std=c99 -O2 _*.c

  -> then run (in a slurm file)
  srun ./ml_convection


HOW TO CREATE A RESTART
  ncks -d time,N out.nc restart.nc

  with N = index of timestep in the file "out.nc"
*/


//#include "grid/multigrid.h"
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/diffusion.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

#define g_ 9.81

/*

TO DO : 
- add dimensions so that they can be saved in the netcdf

*/

/*
DEFAULT PARAMETERS

Dimensions : [Length, Time]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";         // file name of output
// -> Initial conditions
double strat = 0.000002;       // [s-2] N^2 stratification
double Ts = 20.;             // [K] Surface temperature (arbitrary)
// -> Forcing
double qt = 100.;            // [W.m-2] Heat flux
// -> Domain definition
int N_grid = 5 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 2.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;        // dt for output in netcdf
int pad = 4;                // number of 0-padding for ouput files
int nout = 1;               // number of the outfile
char fileout[100];          // name of outfile
// -> physical properties
double nu0 = 0.00025;            // Viscosity for vertical diffusion
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025.;     // [kg.m-3] reference density
double cp = 4.2e3;       // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4;     // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
double D=1.5e-5;
double T0 = 20.;
#define drho(T) (betaT*(T0-T))
#define Tini(z) strat/(g_*betaT)*z + Ts

// double drho0 = 1e-3, T0 = 13.25, T1 = 8., H0 = 100.;
// #define drho(T) (drho0*(T0 - T)/(T0 - T1))
// #define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
//#define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
#include "layered/dr.h"

// diag
double *T_profile;
double dt_mean = 1.;
static FILE * fp;



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
  CFL_H = 1; // Smaller time step


  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  
  //periodic (top);
  //periodic (left);

  // T profile
  T_profile = (double*)calloc(nl, sizeof(double));
  fp  = fopen("T_profile.dat","w"); // reset file
  fclose(fp);

  fprintf (stderr, "Read in parameters!\n");
  run();
  
}

event init(i =  0) {
  geometric_beta (1./3., true); // Varying layer thickness
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
      T[] = Tini(z) -0.01 + (rand() / (RAND_MAX+1.0) * 0.02) ;
      //T[] = T0(z);
      z += h[]/2.;
    }
  }
  create_nc({zb, h, u, w, eta, T}, file_out);
  fprintf (stderr,"Done initialization!\n");

}

// This event adds a heat flux forcing at surface
// event forcing(i++){
//     foreach(){
//     T[0,0,nl-1] = dt*(T[0,0,nl-1] + qt/(rho0*cp*h[0,0,nl-1]));
//   }
// }

// vertical diffusion on T
// This implements a heat flux at the surface
event viscous_term (i++)
{
  foreach()
    vertical_diffusion (point, h, T, dt, D, qt/(D*rho0*cp), 0., 0.);
}

// Writing a 4D netcdf file
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  write_nc();
}

// This event compute layer average of T
event compute_horizontal_avg (t+=dt_mean; t<=tend+1e-10){
  foreach(reduction(+:T_profile[:nl]))
    foreach_layer(){
      T_profile[point.l] += T[] / N;// (N*N); // * dt / dt_mean;
    }
}

// This even writes to a file the layer average
event write_diag(t=0., t+=dt_mean){
    // main worker is writing the file
    if (pid()==0) {
      fp  = fopen("T_profile.dat","a");
      if (fp == NULL){
        fprintf(stderr, "Error opening file T_profile.dat");
        return 2;
      }
      for (int i=0; i<nl; ++i) {
        fprintf (fp, "%f %d %g\n", t, i, T_profile[i]);
      }
      fprintf(fp,"\n");
      fclose(fp);
    }
    // Reset the profile for all workers
    for (int i=0; i<nl; ++i) {
      T_profile[i] = 0.0;
    }
}


// TODO: add a plot of the density (T) profile evolution with time
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
