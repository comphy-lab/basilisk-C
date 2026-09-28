/**
 
 # Wave breaking with stratification (multilayer solver)

*/
const double g_ = 9.81 [1,-2];        // [m.s-2] Gravity

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "hugoj/lib/spectrum.h"
#include "hugoj/lib/diffusionH.h"     // Neumann at top and bot for diffusion
#define PI 3.14159265358979323846

/**
## Default parameters

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char file_out[20] = "out.nc";          // file name of output
// -> Initial conditions
double strat = 0.000002 [0,-2];       // [s-2] N^2 stratification
double Ts = 20. [0,0,1];              // [K] Surface temperature (arbitrary)
double P = 0.02 [1, -1];               // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 5 [];               // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];                   // Number of modes in wavenumber space
int N_power = 5 [];                   // directional spreading coeff
int F_shape = 0 [];                   // shape of the initial spectrum
double kp;         // peak wave number
double Tp;
// -> Forcing
double qt = 100. [-2,-1,0,1];         // [W.m-2] Heat flux at surface
// -> Domain definition
int N_grid = 128;                      // number of x and y gridpoints
int N_layer = 40;                      // number of layers
double L = 200.0 [1];                 // [m] domain size
double h0 = 20.0 [1];                  // [m] depth of water
// -> Runtime parameters
double tend = 30.0 [0,1];              // (x Tp) end time of simulation
// -> saving outputs
double dtout = 5.0 [0,1];             // [s] dt for output in netcdf
double smalltime = 1e-10 [0,1];       // [s] small time increment
// -> physical properties
double Re = 40000;                     // Reynolds number Re = sqrt(g*lambda**3)/nu
double thetaH = 0.503 [0];              // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025. [-3,0,0,0,1];     // [kg.m-3] reference density
double cp = 4.2e3 [0,0,-1,1,-1];      // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4 [0,0,-1];         // [K-1] Thermal expansion coeff for water
double Pr = 1.0;                      // Prandtl number
double diff_T= 1e-5;
double T0 = 20. [0,0,1];              // [°C] Reference temperature

#define drho(T) (betaT*(T0-T))        // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
#define Tini(z) Ts  // Initial temperature profile
#include "layered/dr.h"

int main(int argc, char *argv[])  
{
  kp = 2*PI * coeff_kpL0 / L; // kpL=coeff x 2pi x domain size

  // Settings solver values from namlist values
  L0 = L;
  nu = sqrt(g_*kp)*2*PI/kp/Re; // Re = c*lambda/nu
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  diff_T = nu/Pr;  
  CFL_H = 1; 
  CFL=0.8;
  Tp = 2*PI/sqrt(g_*kp);
  
  tend = tend*Tp;

  // Boundary conditions
  // see also diffusionH.h
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  fprintf (stderr, "Read in parameters!\n");
  run();
  
}

event init(i =  0) {

  /** We read a spectrum using spectrum.h */
  #if WITH_PM
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);  

  /** set eta */
  initial_condition_wave_fft (eta, spectrum, N);
  #endif

  /** set and h*/
  geometric_beta (1., true); // if !=1, varying layer thickness
  foreach(cpu) {
    zb[] = -h0;
    double H = eta[] - zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
    } 
  }

  foreach(cpu) {
    double z = zb[];
    double H = - zb[];
    foreach_layer() {
      z += H*beta[point.l]/2.;
      T[] = Tini(z);
      z += H*beta[point.l]/2.;
    } 
  }

  /** set currents */
  #if WITH_PM
  initial_condition_u_fft (u, spectrum, -h0, 2., N);
  free_spectrum(&spectrum);
  #endif // WITH_PM

  fprintf (stderr,"Done initialization!\n");
    
  create_nc({zb, h, u, w, eta, T}, file_out);
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

/** dump outputs */
event output(t = 0.; t<= tend+smalltime; t+=dtout){
  write_nc();
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
