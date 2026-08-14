#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h"
#define g_ 9.81
#include "hugoj/lib/spectrum.h"

double P = 0.2;               // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 5;               // kpL0 = coeff_kpL0 * pi
int N_mode = 32;                   // Number of modes in wavenumber space
int N_power = 5;                   // directional spreading coeff
int F_shape = 0;                   // shape of the initial spectrum
double kp = PI*10/200.0;         // peak wave number
double Tp;
double h0 = 1.0;                  // [m] depth of water
T_Spectrum spectrum;
char dname[100];

int main(int argc, char *argv[])  
{

  N = argc > 1 ? atoi(argv[1]) : 64;
  nl = argc > 2 ? atoi(argv[2]) : 5;
 
  // Settings solver values from namlist values
  L0 = 200.;
  nu = 0.00025;
  G = g_;
  kp = 2*PI * coeff_kpL0 / L0; // kpL=coeff x 2pi x domain size
  Tp = 2*PI/sqrt(g_*kp);
  
  // Boundary conditions
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  /** initialise once a spectrum */
  spectrum = spectrum_gen_linear(N_mode, N_power, L0, P, kp);
  run();
}

event init_v1(i =  0) {

  /** set eta and h*/
  foreach(cpu) {
    zb[] = -h0;
    eta[] = wave_v1(x, y, spectrum);
    double H = eta[] - zb[];
    foreach_layer() {
      h[] = H/nl;
    } 
  }

    // /** set currents */
    // foreach(cpu) {
    //   double z = zb[];
    //   foreach_layer() {
    //     z += h[]/2.;
    //     coord current = wave_u_v1(x,y,z,spectrum);
    //     u.x[] = current.x; 
    //     u.y[] = current.y; 
    //     w[] = current.z; 
    //     z += h[]/2.;
    //   }
    // }
  sprintf (dname, "out_v1_N%d_nl%d.nc", N, nl);
  create_nc({zb, h, u, w, eta}, dname);
  write_nc();

  fprintf (stderr,"Done initialization v1!\n");
}

event init_v2(i =  0; t<=0.) {

  /** Initialise eta */
  initial_condition_wave_fft (eta, spectrum, N);

  /** set h*/
  foreach(cpu) {
    zb[] = -h0;
    double H = eta[] - zb[];
    foreach_layer() {
      h[] = H/nl;
    } 
  }

  // /** set currents */
  // foreach(cpu) {
  //   double z = zb[];
  //   foreach_layer() {
  //     z += h[]/2.;
  //     coord current = wave_u(x,y,z,spectrum);
  //     u.x[] = current.x; 
  //     u.y[] = current.y; 
  //     w[] = current.z; 
  //     z += h[]/2.;
  //   }
  // }
  
  sprintf (dname, "out_v2_N%d_nl%d.nc", N, nl);
  create_nc({zb, h, u, w, eta}, dname);
  write_nc();

  fprintf (stderr,"Done initialization v2!\n");
  return 1;
}



