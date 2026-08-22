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
#include "hugoj/lib/diffusionH.h"     // Neumann at top and bot for diffusion

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
int coeff_kpL0 = 5 [];               // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];                   // Number of modes in wavenumber space
int N_power = 5 [];                   // directional spreading coeff
int F_shape = 0 [];                   // shape of the initial spectrum
double kp = PI*10/200.0 [-1];         // peak wave number
double Tp;
// -> Forcing
double qt = 100. [-2,-1,0,1];         // [W.m-2] Heat flux at surface
// -> Domain definition
int N_grid = 64;                      // number of x and y gridpoints
int N_layer = 5;                      // number of layers
double L = 200.0 [1];                 // [m] domain size
double h0 = 1.0 [1];                  // [m] depth of water
// -> Runtime parameters
double tend = 2.0 [0,1];              // (x Tp) end time of simulation
double tstart_diag = 1.0;             // (x Tp) start of Stokes drift diag
double tend_diag = 2.0;               // (x Tp) end of Stokes drift diag
double t_iniT = 1.0;                  // (x Tp) time where temperature is initialised
// -> saving outputs
double dtout = 2.0 [0,1];             // [s] dt for output in netcdf
double smalltime = 1e-10 [0,1];       // [s] small time increment
// -> physical properties
double Re = 1000;                     // Reynolds number Re = sqrt(g*lambda**3)/nu
double nu0 = 0.00025 [2,-1];          // Viscosity for vertical diffusion of momentum
double thetaH = 0.5 [0];              // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025. [-3,0,0,0,1];     // [kg.m-3] reference density
double cp = 4.2e3 [0,0,-1,1,-1];      // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4 [0,0,-1];         // [K-1] Thermal expansion coeff for water
double Pr = 1.0;                      // Prandtl number
double diff_T= 1e-5;
double T0 = 20. [0,0,1];              // [°C] Reference temperature
double Trand = 0.1 [0,0,1];           // [°C] Random temperature perturbution

#define drho(T) (betaT*(T0-T))        // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
#define Tini(z) strat/(g_*betaT)*z + Ts  // Initial temperature profile
#include "layered/dr.h"

// diag
double *T_profile;
double *u_profile;
double *zm;
double *hm;
double dt_mean = 1.;

static FILE * fp1;
static FILE * fp2;

int main(int argc, char *argv[])  
{
  /** Building a 'params' array with all parameters from the namlist */
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("tstart_diag", &tstart_diag, "double");
  add_param("tend_diag", &tend_diag, "double");
  add_param("t_iniT", &t_iniT, "double");
  add_param("Re", &Re, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("dtout", &dtout, "double");
  add_param("strat", &strat, "double");
  add_param("Ts", &Ts, "double");
  add_param("rho0", &rho0, "double");
  add_param("cp", &cp, "double");
  add_param("betaT", &betaT, "double");
  add_param("Pr", &Pr, "double");
  add_param("qt", &qt, "double");
  add_param("dt_mean", &dt_mean, "double");

  kp = 2*PI * coeff_kpL0 / L; // kpL=coeff x 2pi x domain size
  
  /** Search for the configuration file with a given path or read params.in */
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);

  // Settings solver values from namlist values
  L0 = L;
  nu0 = sqrt(g_*pow(2*PI/kp, 3))/Re;
  //nu = nu0;
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  diff_T = nu0/Pr;
  CFL_H = 1; 
  CFL=0.8;
  Tp = 2*PI/sqrt(g_*kp);
  
  // Diagnostic times
  tend = tend*Tp;
  tstart_diag = tstart_diag*Tp;
  tend_diag = tend_diag*Tp;
  t_iniT = t_iniT*Tp;

  // Boundary conditions
  // see also diffusionH.h
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  // diags profiles
  T_profile = (double*)calloc(nl, sizeof(double));
  u_profile = (double*)calloc(nl, sizeof(double));
  
  fp1  = fopen("T_profile.dat","w"); // reset file
  fclose(fp1);
  fp2  = fopen("u_profile.dat","w"); // reset file
  fclose(fp2);

  // Stokes drift diag arrays
  zm = (double*)calloc(nl, sizeof(double));
  hm = (double*)calloc(nl, sizeof(double));

  fprintf (stderr, "Read in parameters!\n");
  run();
  
}

event init(i =  0) {

  if (!restore ("restart")) {
    // Warning: no check is done to ensure that you restart with the same
    // resolution/parameters ! Basilisk just reload what is in the dump file.
    
    /** We read a spectrum using spectrum.h */
    T_Spectrum spectrum;
    spectrum = read_spectrum(N_mode);

    /** set eta */
    initial_condition_wave_fft (eta, spectrum, N);

    /** set eta and h*/
    geometric_beta (1/3., true); // if !=0, varying layer thickness
    foreach(cpu) {
      zb[] = -h0;
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
        z += h[]/2.;
        T[] = Tini(z);
        z += h[]/2.;
      } 
    }

    /** set currents */
    initial_condition_u_fft (u, spectrum, -h0, 2., N);

    // initializing diag arrays
    //T_profile[0] = Trand; // <- this passes the dimensional analysis
    //dimensional (T_profile[0] == Trand); // <- this doesnt
    for (int i=0; i<nl; ++i) {
      T_profile[i] = Trand*0.;
      u_profile[i] = L0/DT*0.;
    }
    
    fprintf (stderr,"Done initialization!\n");
    free_spectrum(&spectrum);
  }
  else {
  // Restarting from file 'restart'
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method
      dtmax = 0.01;
      dt = dtnext (dtmax);
    fprintf (stderr, "Restarting ...\n");
  }
  create_nc({zb, h, u, w, eta, T}, file_out);
  write_nc();
 
}

/** vertical diffusion on T and u,w */
event viscous_term (i++; t<=tend)
{
  foreach()
    // vertical_diffusion2 (point, h, T, dt, diff_T, qt/(diff_T*rho0*cp), strat/(g_*betaT));
    vertical_diffusion2 (point, h, T, dt, diff_T, 0., 0.);
  horizontal_diffusion ({T}, diff_T, dt);
  horizontal_diffusion ({u, w}, nu, dt);
}

double* l_avg(scalar var, double* profile){
  /*
  This function computes the layer average of var.

  INPUTS:
    var: scalar (C Basilisk), the variable to average.
    profile: array of double of length nl
  OUTPUTS:
    an array of double with the average values inside
   */
  foreach(reduction(+:profile[:nl]), cpu)
    foreach_layer(){
    #if dimension==1
      profile[point.l] += var[] / N;// (N*N); // * dt / dt_mean;
    #else
      profile[point.l] += var[] / (N*N);
    #endif
    }
  return profile;
}

/** Initialise T after wave spinup */
event initT(t=t_iniT){
  // Compute z mean
  double etam=0.;
  foreach(reduction(+:etam)){
    etam += eta[];
  }
  etam = etam / (N*N);
  hm = l_avg(h, hm); 
  zm[0] = -h0;
  for (int k=1; k < nl; k++){
    zm[k] = zm[k-1] + hm[k-1]/2. + hm[k]/2.;
    fprintf(stderr, "zm[%d]=%g\n", k, zm[k]);
  }
  fprintf(stderr, "<eta>=%g\n",etam);
  // Set T at z mean
  foreach() {
    foreach_layer() {
      T[] = Tini(zm[point.l]); 
    } 
  }
}

#if DUMPS
/** dump outputs */
event output(t = 0.; t<= tend+smalltime; t+=dtout){
  write_nc();
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}
#endif // DUMPS


/** dump for Stokes drift */
#if US_DIAG
event dump_diag_start(t=tstart_diag){
  char dname[100];
  sprintf (dname, "dump_diagstart");
  dump(dname);
}

event output(i++, t<=tend){
  if (t+dt >=tstart_diag && t+dt <= tend_diag){
    write_nc();
  }
}
#endif //US_DIAG

#if DIAG
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





event cleanup(t=end){
  free(zm);
  free(hm);
}

event final_dump(t=end){
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

/**
## TODO:
- 
*/


/**
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
