/**
 
 # Wave breaking with wind forcing (multilayer solver)

*/
const double g_ = 9.81;        // [m.s-2] Gravity


#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "hugoj/lib/spectrum.h"       // Initial conditions generation

/**
## Default parameters

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char file_out[20] = "out.nc"; // file name of output
// -> Initial conditions
double strat = 0.000002;      // [s-2] N^2 stratification
double Ts = 20.;              // [K] Surface temperature (arbitrary)
double P = 0.2;               // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 5;           // kpL0 = coeff_kpL0 * pi
int N_mode = 32;              // Number of modes in wavenumber space
int N_power = 5;              // directional spreading coeff
int F_shape = 0;              // shape of the initial spectrum
double kp;                    // peak wave number
double Tp;
// -> Forcing
double p0 = 1.;               // amplitude of windforcing
// -> Domain definition
int N_grid = 256;              // number of x and y gridpoints
int N_layer = 15;              // number of layers
int N_zlayer = 32;            // for fft init (must be power of 2)
double L = 200.0;             // [m] domain size
double h0 = 100.0;              // [m] depth of water
// -> Runtime parameters
double tend = 100.0;            // (x Tp) end time of simulation
// -> saving outputs
double dtout = 5.0;           // [s] dt for output in netcdf
double smalltime = 1e-10;     // [s] small time increment
// -> physical properties
double Re = 40000;            // Reynolds number Re = sqrt(g*lambda**3)/nu
double thetaH = 0.503;        // theta_h for dumping fast barotropic modes


int main(int argc, char *argv[])  
{
  kp = 2*PI * coeff_kpL0 / L; // kpL=coeff x 2pi x domain size
  
  // Settings solver values from namlist values
  L0 = L;
  nu = sqrt(g_*kp)*2*PI/kp/Re; // Re = cp.lambdap/nu = omegap/kp*lampbdap/nu
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 
  CFL=0.8;
  Tp = 2*PI/sqrt(g_*kp);
  
  // Diagnostic times
  tend = tend*Tp;

  // Boundary conditions
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);
  
  run();
}

event init(i =  0) {
    
  /** We generate a spectrum using spectrum.h */
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);  

  /** set eta */
  initial_condition_wave_fft (eta, spectrum, N);

  /** set h*/
  geometric_beta (1/3., true); // if !=0, varying layer thickness
  foreach(cpu) {
    zb[] = -h0;
    double H = eta[] - zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
    } 
  }


  /** set currents */
  initial_condition_u_fft (u, spectrum, -h0, 3., N, N_zlayer);

  fprintf (stderr,"Done initialization!\n");
  free_spectrum(&spectrum);
  create_nc({zb, h, u, w, eta}, file_out);
  write_nc();
 
}

/** vertical diffusion on u,w */
event viscous_term (i++; t<=tend)
{
  horizontal_diffusion ({u, w}, nu, dt);
}

#if FORCING
#warning "FORCING accel branch active"
event acceleration(i++, last){
  foreach_face(x) {
    double detadx = 0.;
    detadx = (eta[] - eta[-1])/Delta;
    ha.x[0,0,nl-1] += hf.x[]*(p0 * detadx * sin(atan(detadx))); //   
  }
}
#endif

#if DUMPS
/** dump outputs */
event output(t = 0.; t<= tend+smalltime; t+=dtout){
  write_nc();
  // char dname[100];
  // sprintf (dname, "dump_t%g", t);
  // dump(dname);
}
#endif // DUMPS

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
