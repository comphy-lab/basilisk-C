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
double p0 = 50.; // Pa, amplitude of windforcing
double pp0 = 50.; // first guess
double rho = 1000.;
#if HOLD_FORCING || CST_FORCING
//#define a_baro(eta,i) \
  // (1.0 * ( \
  //     -G*(eta[i] - eta[i-1])/Delta \
  //     - (p0/rho) * \
  //       (eta[i+1] - 2.*eta[i] + eta[i-1])/(Delta*Delta) \
  // ))
//#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (p0/rho)*windmask.x[]*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
//#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (p0/rho)*windmask.x[]*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
#define a_baro_p(f,k) (- (p0/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta))
#endif // HOLD_FORCING || CST_FORCING


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
int N_grid = 512;              // number of x and y gridpoints
int N_layer = 15;              // number of layers
int N_zlayer = 32;            // for fft init (must be power of 2)
double L = 200.;             // [m] domain size
double h0 = 100.0;              // [m] depth of water
// -> Runtime parameters
double tend = 5.0;            // (x Tp) end time of simulation
// -> saving outputs
double dtout = 1.0;           // [s] dt for output in netcdf
double smalltime = 1e-10;     // [s] small time increment
// -> physical properties
double Re = 20000.;            // Reynolds number Re = sqrt(g*lambda**3)/nu
double thetaH = 0.503;        // theta_h for dumping fast barotropic modes

double etavar_i = 0.;
double mean_eta_i = 0.;
double Ek_i = 0.;
double Ep_i = 0.;


// Control
double K_1 = 1. [0,0,0];
double relax_dt = 0. [0,1];

double K_2 = 0. [0,0,0];
double T_f = 1.;
double sum_I = 0.;

double K_3 = 0. [0,0,0];
double last_dE=0.;

// Low pass filter
double dE_f, sum_f, tau_f, sum_abs_f;
// p0 clip
double p0_max = 100. [-1,-2,1];

// Wind acts in x only: qcc rewrites .x -> .y inside the
// foreach_dimension()/foreach_face() passes of nh.h/hydro.h, so
// multiplying the wind term by windmask.x[] makes it vanish in the
// y-passes (windmask.y = 0) while the hydrostatic part stays 2D.
vector windmask[];


#define T0  (2.*PI/sqrt(g_*kp))

void energy_stats (double *mean, double *variance, double *pe, double *ke)
{
  double sum = 0.;
  foreach (reduction(+:sum))
    sum += eta[];
  double m = sum/pow (N, dimension);

  double gpe = 0., ke2 = 0.;
  foreach (reduction(+:gpe) reduction(+:ke2)) {
    gpe += sq (eta[] - m)*dv();
    foreach_layer() {
      double norm2 = sq (w[]);
      foreach_dimension()
        norm2 += sq (u.x[]);
      ke2 += norm2*h[]*dv();
    }
  }
  *mean = m;
  *variance = gpe/sq(L0);  // <(eta-m)^2> since sum of dv() = L0^dimension
  *pe = 0.5*rho*G*gpe;
  *ke = 0.5*rho*ke2;
}


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
  CFL=0.5;
  DT = 0.01;
  NITERMIN = 3;
  Tp = 2*PI/sqrt(g_*kp);
  tend = tend*Tp;

  #if CST_FORCING
  p0 *= 10; // wave growth: p0 = ptilde*pref
  #endif

  /** If the dynamic forcing is used, we use PID correction.
    */
  dE_f = 0.;
  sum_f = 0.;
  sum_abs_f = 0.;
  tau_f = T0/2.;  // timescale for low pass
  
  relax_dt = T0; 
  T_f = 1.0; //* T0;
  double Ku = 1./2./relax_dt; // 1./5./relax_dt
  double Tu = 5.*T0;
  double Td = T0/8.;
  
  K_1 = Ku; // P
  K_2 = Ku/Tu; // I K2 = K_1/Ti
  K_3 = 0.*Ku*Td; //2.*Ku*Tu; // D K3 = K_1*Td

  // Boundary conditions
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  dtmax = DT;

  run();
}

/**
 Here a function is defined: it computes the variance and the mean of the
 surface elevation in this 2D model.
 */
// void eta_stats (double *mean, double *variance)
// {
//   double sum = 0., var = 0.;
//   double Nterms = pow (N, dimension);
//
//   foreach (reduction(+:sum))
//     sum += eta[];
//   double m = sum/Nterms;
//
//   foreach (reduction(+:var))
//     var += sq(eta[] - m);
//
//   *mean = m;
//   *variance = var/Nterms;
// }

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
  initial_condition_u_fft (u, spectrum, -h0, N, N_zlayer);

  create_nc({zb, h, u, w, eta}, file_out);
  free_spectrum(&spectrum);

  /** x-only wind direction mask (see a_baro definition) */
  windmask = new vector;
  foreach (cpu) {
    windmask.x[] = 1.;
    windmask.y[] = 0.;
  }

  //initial energy of wavefield
  energy_stats(&mean_eta_i, &etavar_i, &Ep_i, &Ek_i);

  fprintf(stderr, "Initial Hskp = %f\n", 4*sqrt(etavar_i)*kp);
  fprintf(stderr, "target variance = %f\n", etavar_i);
  fprintf(stderr, "target Hs = %f\n", 4*sqrt(etavar_i));
  fprintf(stderr, "target Et = %g\n", Ep_i+Ek_i);

  fprintf (stderr,"Done initialization!\n");
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
  /**
  IN 2D (x,z):
  input_power = -sum_xy(dp/dx Q)dx
             with p = p0.deta/dx
                  Q = sum_layer(u*h)
      
        so input_power = p0* sum_xy(sum_layer(d2eta/dx2*u*h*dx))


   Energy deficit = rho*g*(variance_eta_target - variance_eta_now)*L0**2

   timescale of relaxation = relax_dt

   we match input power = deficit of relax_dt:

     -sum_xy(dp/dx Q)dx = rho*g*(variance_eta_target - variance_eta_now)*L0**2/relax_dt
  */

  /** Computing the error on wave energy */
  // double etam = 0.;
  // double etavar_current = 0.;
  // eta_stats (&etam, &etavar_current);
  // double dE = rho*G*(etavar_i - etavar_current)*sq(L0);

  double etam = 0., etavar_current = 0., pe = 0., ke = 0.;
  energy_stats (&etam, &etavar_current, &pe, &ke);
  double dE = (Ek_i+Ep_i) - (pe + ke);

  /** Computing input power */
  // P = p0*sum
  double sum = 0.; 
  foreach(reduction(+:sum)) {
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer(){
      // fprintf(stderr, "i %d, j %d, k %d, h %f, u.x %f, etaxx %f\n",
      //     point.i,point.j,point.l,h[],u.x[],etaxx);
      sum += - h[] * u.x[] * etaxx * dv(); 
    }
  }
  /**
  Low pass filter on error and sum
  */
  if (i>0){
    dE_f += (dE  - dE_f ) * dt / tau_f;
    sum_f += (sum - sum_f) * dt / tau_f;
  }
  /**
  safeguard: never divide by near zero sum  
  */
  sum_abs_f += (fabs(sum) - sum_abs_f) * dt / tau_f;
  double sum_den = (sum_f >= 0. ? 1. : -1.) * fmax (fabs(sum_f), 0.1*sum_abs_f);
  /**
  ## PID Control

  p0 = 1/sum * K1( dE + 1/Ti*int_t dE dt + Td*dE/dt )

  */
  double rhs = 0.;
  /**
  P: proportional constrain
  */
  rhs += K_1 * dE_f;

  /**
  I: integral constrain   

  */
  //sum_I += dE_f*dt;
  if (fabs((rhs + K_2/T_f*sum_I)/sum_den) <= p0_max)
    sum_I += dE_f * dt;
  rhs += K_2 / T_f * sum_I;

  /**
  D: derivative constrain
  */  
  double dedt = 0.;  
  dedt = (dE_f-last_dE)/dt;
  rhs += K_3*dedt;
  last_dE = dE_f;
  /**
  Finally we update p0
  */
  double dp0 = rhs/sum_den;
  p0 = pp0 + dp0;

  /**
  p0 limiter
  */
  if (p0 >  p0_max) p0 =  p0_max;
  if (p0 < -p0_max) p0 = -p0_max;

  //p0 = 1.0;

  //p0 = -rho*G*dE*pow(L0,dimension)/relax_dt / sum;
  // fprintf(stderr,
  //       "i=%d dE=%g predicted=%g sum=%g relax_dt=%g p0=%g\n",
  //       i, dE,-p0*sum*relax_dt, sum, relax_dt, p0);
  // if (i%10==0){
  //   fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g sum=%g sum_den=%g p0=%g\n",
  //       i, etam, etavar_current, dE, sum, sum_den, p0);
  // }
  fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g Et=%g sum=%g sum_den=%g p0=%g\n",
        i, etam, etavar_current, dE, pe+ke,  sum, sum_den, p0);

}
#endif  // HOLD_FORCING


// Note
// For now I have to use this face_fields event to add the forcing.
// When i'll be able to use a vector inside the a_baro macro I will 
// be able to use the original machinery with the a_baro macro overloaded.
#if HOLD_FORCING || CST_FORCING
event face_fields (i++,first) {
  foreach_face (x) {
    double ax = a_baro_p (eta, 0);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
    }
  }
}
#endif // HOLD_FORCING || CST_FORCING



/** Log energy evolution */
event logfile (i++; t <= tend) // target: at least 300*T0 
{ 
  if (i%10==0){
    double etam = 0., var_eta=0., ke=0., gpe=0.;
    energy_stats(&etam, &var_eta, &gpe, &ke);
    fprintf (stdout, "%g %g %g\n", t/T0, ke, gpe);
  }

}

/** dump outputs */
event output(i++)
{
  if (i%10==0){
    write_nc();
  }
}
