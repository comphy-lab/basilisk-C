/**
 
 # Wave breaking with wind forcing (multilayer solver)

compute upper incomplete gamma funciton : https://www.azcalculator.com/calculators/incomplete-gamma-function
how to do it in python ? In C ?


generate movie with
visu_3Dmovie.py --speed-factor 5 --fps 1 out.nc movie
*/
const double g_ = 9.81;        // [m.s-2] Gravity

#include "grid/multigrid.h"

// -> Forcing
double p0ref = 0.00625;               // amplitude of windforcing
double p0 = 0.00625; // Pa
double rho = 1000.;
// #if HOLD_FORCING || CST_FORCING
// #define wind_pressure(eta,i)  (p0/rho)*(eta[i+1] - eta[i-1])/(2*Delta)
// #define p_baro(eta,i) (-G*eta[i]) - wind_pressure(eta,i))
// #define a_baro(eta,i) (gmetric(i)*(p_baro(eta,i)-p_baro(eta,i-1))/Delta)

#define a_baro(eta,i) \
  (gmetric(i) * ( \
      -G*(eta[i] - eta[i-1])/Delta \
      - (p0/rho) * \
        (eta[i+1] - 2.*eta[i] + eta[i-1])/(Delta*Delta) \
  ))

// #endif // HOLD_FORCING


#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "hugoj/lib/spectrum.h"       // Initial conditions generation

/**
## Default parameters

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char file_out[20] = "out.nc"; // file name of output
// -> Initial conditions
double P = 0.02;               // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 5;           // kpL0 = coeff_kpL0 * pi
int N_mode = 32;              // Number of modes in wavenumber space
int N_power = 5;              // directional spreading coeff
int F_shape = 0;              // shape of the initial spectrum
double kp=2*PI;                    // peak wave number
double Tp;

// -> Domain definition
int N_grid = 64;              // number of x and y gridpoints
int N_layer = 15;              // number of layers
int N_zlayer = 32;            // for fft init (must be power of 2)
double L = 200.;             // [m] domain size
double h0 = 100.0;              // [m] depth of water
// -> Runtime parameters
double tend = 10.0;            // (x Tp) end time of simulation
// -> saving outputs
double dtout = 5.0;           // [s] dt for output in netcdf
double smalltime = 1e-10;     // [s] small time increment
// -> physical properties
double Re = 40000.;            // Reynolds number Re = sqrt(g*lambda**3)/nu
double thetaH = 0.503;        // theta_h for dumping fast barotropic modes

double relax_dt;
double etavar_i = 0.;
double etavar_current = 0.;

#define T0  (2.*PI/sqrt(g_*kp))

int main(int argc, char *argv[])  
{
  kp = 2*PI * coeff_kpL0 / L; // kpL=coeff x 2pi x domain size

  L0 = L;
  nu = sqrt(g_*kp)*2*PI/kp/Re; // Re = cp.lambdap/nu = omegap/kp*lampbdap/nu
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 
  CFL=0.8;
  Tp = 2*PI/sqrt(g_*kp);
  p0 *= p0ref;
  tend = tend*Tp;

  #if CST_FORCING
  p0 *= 10; // wave growth: p0 = ptilde*pref
  #endif
  relax_dt = 5*T0; 

  // Boundary conditions
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);
  
  run();
}

/**
 Here a function is defined: it computes the variance and the mean of the
 surface elevation in this 2D model.
 */
void eta_stats2D (double *mean, double *variance)
{
  double sum = 0., var = 0.;
  foreach (reduction(+:sum))
    sum += eta[]; 
  *mean = sum/(N*N);
  foreach (reduction(+:var))
    var += sq(eta[] - *mean);
  *variance = var/(N*N); 
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
  double mean_eta_i = 0.;
  //initial energy of wavefield
  eta_stats2D(&mean_eta_i, &etavar_i);
  fprintf(stderr, "Initial Hskp = %f\n", 4*sqrt(etavar_i)*kp);
  fprintf(stderr, "target variance = %f\n", etavar_i);
  fprintf(stderr, "target Hs = %f\n", 4*sqrt(etavar_i));
  

  /** set currents */
  initial_condition_u_fft (u, spectrum, -h0, N, N_zlayer);

  fprintf (stderr,"Done initialization!\n");
  free_spectrum(&spectrum);
  create_nc({zb, h, u, w, eta}, file_out);
}

/** vertical diffusion on u,w */
event viscous_term (i++; t<=tend) { 
  // Vertical diffusion is in diffusion.h and nh.h
  horizontal_diffusion ({u, w}, nu, dt);
}


/**
 ## Wind forcing
*/
#if HOLD_FORCING
event update_p0 (i++)
{
  double etam;
  eta_stats2D (&etam, &etavar_current);
  double dE = (etavar_i - etavar_current);
  double sum = 0.;
  foreach(reduction(+:sum)) {
    double integrated_transport=0.;
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer()
      integrated_transport += h[] * u.x[]; // h * detadx * u
    sum += integrated_transport*etaxx*dv();
  }
  p0 = -G*rho*dE/(relax_dt*sum);
  if (i%10==0){
    fprintf(stderr,
        "i=%d, t=%g, etavar= %g, dE=%g sum=%g dt=%g p0=%g predicted=%g\n",
        i, t, etavar_current, dE, sum, relax_dt, p0,
        -p0*sum*relax_dt/rho);
  }
}
#endif  // HOLD_FORCING


/** dump outputs */
event output(t = 0.; t<= tend+smalltime; t+=dtout)
{
  write_nc();
}

// #if _GPU && SHOW && !_CUDA
// event display (i++){
//   vector u = lookup_vector ("u14");
//   output_ppm (u.x, min = -1.5, max = 1.5, fps = 30, fp = NULL, map=gray);
// }
// #endif

// event final_dump(t=end){
//   char dname[100];
//   sprintf (dname, "dump_t%g", t);
//   dump(dname);
// }



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
