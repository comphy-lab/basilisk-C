/**
 
 # Wave breaking with stratification (multilayer solver)

*/
const double g_ = 9.81 [1,-2];        // [m.s-2] Gravity

//#include "grid/multigrid.h"
//#include "grid/multigrid1D.h"
//#include "grid/cuda/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "hugoj/lib/spectrum.h"       // Initial conditions generation
#include "hugoj/lib/diffusionH.h"     // Neumann at top and bot

/**
## Default parameters

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";          // file name of output
// -> Initial conditions
double strat = 0.000002 [0,-2];       // [s-2] N^2 stratification
double Ts = 20. [0,0,1];              // [K] Surface temperature (arbitrary)
double P = 0.2 [1, -1];               // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];               // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];                   // Number of modes in wavenumber space
int N_power = 5 [];                   // directional spreading coeff
int F_shape = 0 [];                   // shape of the initial spectrum
double kp = PI*10/200.0 [-1];         // peak wave number
double Tp;
// -> Forcing
double qt = 100. [-2,-1,0,1];         // [W.m-2] Heat flux
// -> Domain definition
int N_grid = 5;                       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];                 // domain size
int N_layer = 5;                      // number of layers
double h0 = 1.0 [1];                  // depth of water
// -> Runtime parameters
double tend = 2.0 [0,1];              // end time of simulation
double tstart_diag = 1.0;
double tend_diag = 2.0;
// -> saving outputs
double dtout = 2.0 [0,1];             // dt for output in netcdf
double smalltime = 1e-10 [0,1];             // (s) small time increment
int pad = 4 [0];                      // number of 0-padding for ouput files
int nout = 1 [0];                     // number of the outfile
char fileout[100];                    // name of outfile
// -> physical properties
double Re = 1000;                     // Reynolds number Re = sqrt(g*lambda**3)/nu
double nu0 = 0.00025 [2,-1];          // Viscosity for vertical diffusion
double thetaH = 0.5 [0];              // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025. [-3,0,0,0,1];     // [kg.m-3] reference density
double cp = 4.2e3 [0,0,-1,1,-1];      // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4 [0,0,-1];         // [K-1] Thermal expansion coeff for water
double Diff_T = 1.5e-5 [2,-1];         // [m2.s-1] Scalar vertical diffusion coeff
double T0 = 20. [0,0,1];              // [°C] Reference temperature
double Trand = 0.1 [0,0,1];           // [°C] Random temperature perturbution

#define drho(T) (betaT*(T0-T))        // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
#define Tini(z) strat/(g_*betaT)*z + Ts
#include "layered/dr.h"

// diag
double *T_profile;
double *u_profile;
double dt_mean = 1.;

static FILE * fp1;
static FILE * fp2;
scalar T_ini[];


int main(int argc, char *argv[])  
{
 
  /** Building a 'params' array with all parameters from the namlist */
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("Re", &Re, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("dtout", &dtout, "double");
  add_param("strat", &strat, "double");
  add_param("Ts", &Ts, "double");
  add_param("rho0", &rho0, "double");
  add_param("cp", &cp, "double");
  add_param("betaT", &betaT, "double");
  add_param("qt", &qt, "double");
  add_param("dt_mean", &dt_mean, "double");

  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
  /** Search for the configuration file with a given path or read params.in */
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);

  // Settings solver values from namlist values
  L0 = L;
  nu0 = sqrt(g_*pow(2*PI/kp, 3))/Re;
  nu = nu0;
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 
  CFL=0.8;
  Tp = 2*PI/sqrt(g_*kp);
  
  tend = 13*Tp;
  tstart_diag = 10*Tp;
  tend_diag = 12*Tp;

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  // diags
  T_profile = (double*)calloc(nl, sizeof(double));
  u_profile = (double*)calloc(nl, sizeof(double));
  fp1  = fopen("T_profile.dat","w"); // reset file
  fclose(fp1);
  fp2  = fopen("u_profile.dat","w"); // reset file
  fclose(fp2);

  fprintf (stderr, "Read in parameters!\n");
  

  run();
  
}

event init(i =  0) {
  T_ini = new scalar[nl];
  geometric_beta (1/3., true); // if !=0, varying layer thickness
  
  /** We read a spectrum using spectrum.h */
  T_Spectrum spectrum;
  spectrum = read_spectrum(N_mode);

  /** set eta and h*/
  foreach(cpu) {
    zb[] = -h0;
    eta[] = wave(x, y, spectrum);
    double H = eta[] - zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
    } 
  }

  /** set a temporary T_ini field that will be used to initialse T at a later
    time */
  foreach(cpu) {
    double z = zb[];
    foreach_layer() {
      z+=h0*beta[point.l]/2.;
      T_ini[] = Tini(z); // + noise();
      T[] = Tini(z);
      z+=h0*beta[point.l]/2.;
    } 
  }

  /** set currents */
  foreach(cpu) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      coord current = wave_u(x,y,z,spectrum);
      u.x[] = current.x; 
      u.y[] = current.y; 
      w[] = current.z; 
      z += h[]/2.;
    }
  }

  // initializing diag arrays
  //T_profile[0] = Trand; // <- this passes the dimensional analysis
  //dimensional (T_profile[0] == Trand); // <- this doesnt
  for (int i=0; i<nl; ++i) {
    T_profile[i] = Trand*0.;
    u_profile[i] = L0/DT*0.;
  }
  create_nc({zb, h, u, w, eta, T}, file_out);
  fprintf (stderr,"Done initialization!\n");
}


/** Initialise T after wave spinup 
 
This is a preliminary way of doing the T initialisation.
I see that on the T field (after about 15s=3*Tp), the left boundary does some
strange things. tested at N=128 nl=30. I'm testing on N=512 to see if it is
still the case.

Maybe a problem with boundary condition because I initialise after the fact T ?
Do I need to reset T to 0 at i=1 or i=2, then initialise at T=3Tp ?

 
  */
event initT(t=3*Tp){
  foreach() {
    foreach_layer() {
      T[] = T_ini[]; // + noise();
    } 
  }
}

/** dump outputs */
// event output(t = 0.; t<= tend+smalltime; t+=dtout){
//   write_nc();
//   char dname[100];
//   sprintf (dname, "dump_t%g", t);
//   dump(dname);
// }


/** dump for Stokes drift */
#if US_DIAG
event output(i++, t<=tend){
  if (t+dt >=0. && t+dt <= Tp){
    write_nc();
  }
}
#endif US_DIAG


#if DIAG
double* l_avg(scalar var, double* profile){
  /*
  This function computes the layer average of var.

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
event compute_layer_avg (t+=dt_mean; t<=tend+smalltime){

  T_profile = l_avg(T, T_profile);
  u_profile = l_avg(u.x, u_profile);
}



// This even writes to a file the layer average
event write_diag(t=0.; t+=dt_mean){

    write_profile("T_profile.dat", T_profile, fp1);
    write_profile("u_profile.dat", u_profile, fp2);

    // Reset the profile for all workers
    for (int i=0; i<nl; ++i) {
      T_profile[i] = 0.0;
      u_profile[i] = 0.0;
    }
}


event cleanup(t=end){
  free(T_profile);
  free(u_profile);
}
#endif // DIAG

event final_dump(t=end){
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

/**
## TODO:

- impose initial stratification at bottom for T? -> if I impose a flux at top
yes, else no. 
*/


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

~~~pythonplot Temperature profile evolution
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("T_profile.dat")
nl=30
nt = data.shape[0]//nl

fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles.png", dpi=150)
plt.show()
~~~

**/
